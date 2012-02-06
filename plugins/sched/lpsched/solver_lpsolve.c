/*****************************************************************************\
 *  solver_lpsolve.c - linear programming based scheduler plugin.
 *
 *  Developed by Seren Soner and Can Ozturan 
 *  from Bogazici University, Istanbul, Turkey
 *  as a part of PRACE WP9.2
 *  Contact seren.soner@gmail.com, ozturaca@boun.edu.tr if necessary
 
 *  Used backfill plugin as base code. There are some functions that are
    left untouched.

 *  For the theory behind the plugin, see the manuscript
 *  doi: ....
\*****************************************************************************/

#include "solver_lpsolve.h"

extern int solve_allocation(int m, int n, int timeout,
		sched_nodeinfo_t *node_array, solver_job_list_t *job_array)
{
	lprec *lp;
	solver_job_list_t *sjob_ptr;
	struct job_details *job_det_ptr;
	int Ncol = (2 * m + 2) * n;
	int i, j, k;
	int *colno = (int *)malloc(2 * n * sizeof(*colno));
	REAL *row = (REAL *)malloc(2 * n * sizeof(*row));
	REAL *var = (REAL *)malloc(Ncol * sizeof(*var));
	REAL *sparserow;
	int *sparsecol;
	char name[20];
	/*
	m = nodesize
	n = windowSize
	s_(n) binary variable, selected or not : 1 - n
	x_(m*n) on_node cpu assignment variables : n + 1 : n*(m+1)
	t_(n*m) whether a job uses a node or not : n*(m+1)+1 : n*(2m+1)
	c_(n) cost variable : (n*2m+1) + 1 : n*(2m+2)
	sum: n*(2m+2)
	*/
	lp = make_lp(0,Ncol);
	set_timeout(lp, timeout);
	if (lp == NULL)
		return -1;

	for (j = 0; j < n; j++) {
		sprintf(name, "s_%d",j+1);
		set_col_name(lp, j + 1, name);
	}
	for (i = 0; i < m; i++) {
		for (j = 0; j < n; j++) {		
			sprintf(name, "x_%d_%d",i+1,j+1);
			set_col_name(lp, (1 + i) * n + 1 + j, name);
		}
	}
	for (j = 0; j < n; j++) {		
		sprintf(name, "c_%d",j+1);
		set_col_name(lp, n * (2 * m + 1) + j + 1, name);
		for (i = 0; i < m; i++) {
			sprintf(name, "t_%d_%d",j+1,i+1);
			set_col_name(lp, n * (m + 1) + j * m + 1 + i, name);
		}
	}

	set_add_rowmode(lp, TRUE);
        /*
	sum over nodes (allocated cpu) should be equal to job's requested cpu
	sum_i(x_ij) == r_j * s_j
	*/
	sparserow = (REAL*)malloc((m + 1) * sizeof(*sparserow));
	sparsecol = (int*)malloc((m + 1) * sizeof(*colno));
	for (j = 0; j < n; j++) {
		for (i = 0; i < m; i++) {
			sparserow[i] = 1.0;
			sparsecol[i] = (1 + i) * n + j + 1;
		}
		sparserow[m] = (int)(-job_array[j].min_cpus);
		sparsecol[m] = j + 1;
		add_constraintex(lp, m + 1, sparserow, sparsecol, EQ, 0);
	}
	free(sparserow);
	free(sparsecol);

        /*
	sum over jobs for cpu should be available 
	sum_j(x_ij) <= R_i
	*/
	sparserow = (REAL*)malloc(n * sizeof(*sparserow));
	sparsecol = (int*)malloc(n * sizeof(*colno));
	for (i = 0; i < m; i++) {
		for (j = 0; j < n; j++) {
			sparserow[j] = 1.0;
			sparsecol[j] = (1 + i) * n + j + 1;
		}
		add_constraintex(lp, n, sparserow, sparsecol, LE, 
			(int)(node_array[i].rem_cpus));
	}
	free(sparserow);
	free(sparsecol);

        /*
	sum over jobs for gpu should be available on that node
	sum_j(t_ji * g_j) <= G_i
	*/
	sparserow = (REAL*)malloc(n * sizeof(*sparserow));
	sparsecol = (int*)malloc(n * sizeof(*colno));
	for (i = 0; i < m; i++) {
		for (j = 0; j < n; j++) {
			sparserow[j] = job_array[j].gpu;
			sparsecol[j] = n * (m + 1) + 1 + i * n + j;
		}
		add_constraintex(lp, n, sparserow, sparsecol, LE, 
			(int)(node_array[i].rem_gpus));
	}
	free(sparserow);
	free(sparsecol);

	/* 
	FUTURE WORK FOR TOPOLOGY AWARENESS
	t_ji is 1 if job j is allocated resourced in node i
	t_ji = 1 if x_ij > 0
	*/
	sparserow = (REAL*)malloc(2 * sizeof(*sparserow));
	sparsecol = (int*)malloc(2 * sizeof(*colno));
	for (j = 0; j < n; j++) {
		for (i = 0; i < m; i++) {
			sparserow[0] = 1.0;
			sparsecol[0] = (1 + i) * n + 1 + j;
			sparsecol[1] = n * (m + 1) + 1 + j * m + i;
			sparserow[1] = -1.0*MIN((int)(job_array[j].min_cpus),
				(int)(node_array[i].rem_cpus));
			add_constraintex(lp, 2, sparserow, sparsecol, LE, 0);
			sparserow[0] = 1.0;
			sparsecol[0] = (1 + i) * n + 1 + j;
			sparsecol[1] = n * (m + 1) + 1 + j * m + i;
			sparserow[1] = -1;
			add_constraintex(lp, 2, sparserow, sparsecol, GE, 0);
		}
	}
	free(sparserow);
	free(sparsecol);

        /*
	cost is only selected nodes 
	c_j constraint; sum_i(t_ji) / (m + 1)
	*/
	sparserow = (REAL*)malloc((m + 1) * sizeof(*sparserow));
	sparsecol = (int*)malloc((m + 1) * sizeof(*colno));
	for (j = 0; j < n; j++) {
		for (i = 0; i < m; i++) {
			sparserow[i] = 1.0;
			sparsecol[i] = n * (m + 1) + 1 + j * m + i;
		}
		sparserow[m] = -(m + 1);
		sparsecol[m] = n * (2 * m + 1) + j + 1;
		add_constraintex(lp, m + 1, sparserow, sparsecol, EQ, 0);
	}
	free(sparserow);
	free(sparsecol);

        /*
	min_nodes <= c_j * (m + 1) <= max_nodes
	*/
	sparserow = (REAL*)malloc((2) * sizeof(*sparserow));
	sparsecol = (int*)malloc((2) * sizeof(*colno));
	for (j = 0; j < n; j++) {
		sparserow[0] = 1.0 + 1.0 * m;
		sparserow[1] = -(int)(job_array[j].min_nodes);
		sparsecol[0] = n * (2 * m + 1) + 1 + j;
		sparsecol[1] = j + 1;
		add_constraintex(lp, 2, sparserow, sparsecol, GE, 0);
		sparserow[0] = 1.0 + 1.0 * m;
		sparserow[1] = -(int)(job_array[j].max_nodes);
		sparsecol[0] = n * (2 * m + 1) + 1 + j;
		sparsecol[1] = j + 1;
		add_constraintex(lp, 2, sparserow, sparsecol, LE, 0);
	}
	free(sparserow);
	free(sparsecol);

	set_add_rowmode(lp, FALSE);
	/* set the objective function to p_j * (s_j - c_j) */
	for (j = 0; j < n; j++) {
		colno[j] = j + 1;
		row[j] = (double)(job_array[j].priority);
		colno[j + n] = j + 1 + n;
		row[j + n] = -job_array[j].priority;		
	}
	set_obj_fnex(lp, n, row, colno);
	set_maxim(lp);
	/* 
	set the s_j to be binary
	x_ij to integer
	t_ij to binary
	*/	
	for (j = 1; j <= n; j++) {
		set_int(lp, j, TRUE); /* s_j */
		set_bounds(lp, j, 0.0, 1.0);
	}
	for (j = n + 1; j <= n * (m + 1); j++) {
		set_int(lp, j, TRUE); /* x_ij */
		set_bounds(lp, j, 0.0, get_infinite(lp));
	}
	for (j = n * (m + 1) + 1; j <= n * (2 * m + 1); j++) {
		set_int(lp, j, TRUE); /* t_ij */
		set_bounds(lp, j, 0.0, 1.0);
	}
	set_verbose(lp,NORMAL);
	/*	
	write_lp(lp, "model.lp");
	*/
	k = solve(lp);
	get_variables(lp, var);	
	delete_lp(lp);
	for (j = 0; j < n; j++) {
		if (var[j] > 0) {
			sjob_ptr = &job_array[j];
			job_det_ptr = sjob_ptr->job_ptr->details;
			sjob_ptr->node_bitmap = (bitstr_t *) 
				bit_alloc (node_record_count);
			job_det_ptr->req_node_bitmap = (bitstr_t *) 
				bit_alloc (node_record_count);
			sjob_ptr->onnodes = (uint32_t *) xmalloc 
				(sizeof(uint32_t) * node_record_count);
			job_det_ptr->req_node_layout = (uint16_t *) xmalloc
				(sizeof(uint16_t) * node_record_count);
			job_det_ptr->req_node_bitmap = (bitstr_t *) bit_alloc
				(node_record_count);
			for (i = 0; i < m; i++) {
				k = (1 + i) * n + j;
				if (var[k] > 0) {
					bit_set (sjob_ptr->node_bitmap, 
						(bitoff_t) (i));
					bit_set (job_det_ptr->req_node_bitmap,
						(bitoff_t) (i));		
					node_array[i].rem_cpus -= var[k];
					node_array[i].rem_gpus -= sjob_ptr->gpu;
					sjob_ptr->onnodes[i] = var[k];
					job_det_ptr->req_node_layout[i] = 
						sjob_ptr->onnodes[i]; 
					sjob_ptr->alloc_total += var[k];
				}
			}
		} else
			job_array[j].alloc_total = 0;
	} 
	return k;
}

