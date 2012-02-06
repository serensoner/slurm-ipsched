#include <stdlib.h>
#include <stdio.h>
#include "lp_lib.h"
#include "lpsched.h"

/* 
returns -1 if cannot construct a model
	0, 1, 2... if solve_lp outputs
*/
int solvewindow(int nodeSize, int windowSize, int timeout,
		sched_nodeinfo_t *node_array, solver_job_list_t *job_array)
{
	lprec *lp;
	int Ncol = (2 * nodeSize + 2) * windowSize + 1;
	int i, j, k;
	int *colno = (int *) malloc(2 * windowSize * sizeof(*colno));
	REAL *row = (REAL *)malloc(2 * windowSize * sizeof(*row));
	REAL *var = (REAL *)malloc(Ncol * sizeof(*var));
	REAL *sparserow;
	int *sparsecol;
	char varname[20];
	/*
	m = nodesize
	n = windowsize
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

	for (j = 0; j < windowSize; j++) {
		sprintf(varname, "s_%d",j+1);
		set_col_name(lp, j + 1, varname);
	}
	for (i = 0; i < nodeSize; i++) {
		for (j = 0; j < windowSize; j++) {		
			sprintf(varname, "x_%d_%d",i+1,j+1);
			set_col_name(lp, (1 + i) * windowSize + 1 + j, varname);
		}
	}
	for (j = 0; j < windowSize; j++) {		
		sprintf(varname, "c_%d",j+1);
		set_col_name(lp, windowSize * (2 * nodeSize + 1) + j + 1, varname);
		for (i = 0; i < nodeSize; i++) {
			sprintf(varname, "t_%d_%d",j+1,i+1);
			set_col_name(lp, windowSize * (nodeSize + 1) + j * nodeSize + 1 + i, varname);
		}
	}

	set_add_rowmode(lp, TRUE);
        /*
	sum over nodes (allocated cpu) should be equal to job's requested cpu
	sum_i(x_ij) == r_j * s_j
	*/
	sparserow = (REAL*)malloc((nodeSize + 1) * sizeof(*sparserow));
	sparsecol = (int*)malloc((nodeSize + 1) * sizeof(*colno));
	for (j = 0; j < windowSize; j++) {
		for (i = 0; i < nodeSize; i++) {
			sparserow[i] = 1.0;
			sparsecol[i] = (1 + i) * windowSize + j + 1;
		}
		sparserow[nodeSize] = -job_array[j].min_cpus;
		sparsecol[nodeSize] = j + 1;
		add_constraintex(lp, nodeSize + 1, sparserow, sparsecol, EQ, 0);
	}
	free(sparserow);
	free(sparsecol);

        /*
	sum over jobs for cpu should be available 
	sum_j(x_ij) <= R_i
	*/
	sparserow = (REAL*)malloc(windowSize * sizeof(*sparserow));
	sparsecol = (int*)malloc(windowSize * sizeof(*colno));
	for (i = 0; i < nodeSize; i++) {
		for (j = 0; j < windowSize; j++) {
			sparserow[j] = 1.0;
			sparsecol[j] = (1 + i) * windowSize + j + 1;
		}
		add_constraintex(lp, windowSize, sparserow, sparsecol, LE, node_array[i].rem_cpus);
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
	for (j = 0; j < windowSize; j++) {
		for (i = 0; i < nodeSize; i++) {
			sparserow[0] = 1.0;
			sparsecol[0] = (1 + i) * windowSize + 1 + j;
			sparsecol[1] = windowSize * (nodeSize + 1) + 1 + j * nodeSize + i;
			sparserow[1] = -job_array[j].min_cpus;
			add_constraintex(lp, 2, sparserow, sparsecol, LE, 0);
			sparserow[1] = -node_array[i].rem_cpus;
			add_constraintex(lp, 2, sparserow, sparsecol, LE, 0);
		}
	}
	free(sparserow);
	free(sparsecol);
        /*
	cost is only selected nodes 
	c_j constraint; sum_i(t_ji) / (nodeSize + 1)
	*/
	sparserow = (REAL*)malloc((nodeSize + 1) * sizeof(*sparserow));
	sparsecol = (int*)malloc((nodeSize + 1) * sizeof(*colno));
	for (j = 0; j < windowSize; j++) {
		for (i = 0; i < nodeSize; i++) {
			sparserow[i] = 1.0;
			sparsecol[i] = windowSize * (nodeSize + 1) + 1 + j * nodeSize + i;
		}
		sparserow[nodeSize] = -(nodeSize + 1);
		sparsecol[nodeSize] = windowSize * (2 * nodeSize + 1) + j + 1;
		add_constraintex(lp, nodeSize + 1, sparserow, sparsecol, EQ, 0);
	}
	free(sparserow);
	free(sparsecol);

        /*
	min_nodes <= c_j * (nodeSize + 1) <= max_nodes
	*/
	sparserow = (REAL*)malloc((1) * sizeof(*sparserow));
	sparsecol = (int*)malloc((1) * sizeof(*colno));
	for (j = 0; j < windowSize; j++) {
		sparserow[0] = 1.0 + 1.0 * nodeSize;
		sparsecol[0] = windowSize * (2 * nodeSize + 1) + 1 + j;
		add_constraintex(lp, 1, sparserow, sparsecol, GE, job_array[j].min_nodes);
		add_constraintex(lp, 1, sparserow, sparsecol, LE, job_array[j].max_nodes);
	}
	free(sparserow);
	free(sparsecol);

	set_add_rowmode(lp, FALSE);
	/* set the objective function to p_j * (s_j - c_j) */
	for (j = 0; j < windowSize; j++) {
		colno[j] = j + 1;
		row[j] = job_array[j].priority;
		/* topo-aware 
		colno[j + windowSize] = j + 1 + windowSize;
		row[j + windowSize] = -job_array[j].priority;		
		*/
	}
	set_obj_fnex(lp, windowSize, row, colno);
	set_maxim(lp);
	/* 
	set the s_j to be binary
	x_ij to integer
	t_ij to binary
	*/	
	for (j = 1; j <= windowSize; j++) {
		set_int(lp, j, TRUE); /* s_j */
		set_bounds(lp, j, 0.0, 1.0);
	}
	for (j = windowSize + 1; j <= windowSize * (nodeSize + 1); j++) {
		set_int(lp, j, TRUE); /* x_ij */
		set_bounds(lp, j, 0.0, get_infinite(lp));
	}
	for (j = windowSize * (nodeSize + 1) + 1; j <= windowSize * (2 * nodeSize + 1); j++) {
		set_int(lp, j, TRUE); /* t_ij */
		set_bounds(lp, j, 0.0, 1.0);
	}
	set_verbose(lp,NORMAL);
	write_lp(lp, "model.lp");
	set_scaling(lp,SCALE_MEAN);
	k = solve(lp);
	print_solution(lp, 1);
	get_variables(lp, var);	
	for (j = 0; j < windowSize; j++) {
		if (var[j] > 0) {
			printf("%d: ", j + 1);
			for (i = 0; i < nodeSize; i++) {
				k = (1 + i) * windowSize + j;
				if (var[k] > 0)
					printf("%d/%d ",(int)var[k],i + 1);
			}
			printf("\n");
		}
	} 
	return k;
}

int main() 
{
	int m = 8;
	int n = 3;
	int i, j;
	/* need improvement here */
	/* first see if node has rem_cpu > 0 */
	/* if not, do not include in node_array[] */
	sched_nodeinfo_t *node_array = (sched_nodeinfo_t*)malloc(m * sizeof(sched_nodeinfo_t));
	solver_job_list_t *job_array = (solver_job_list_t*)malloc(n * sizeof(solver_job_list_t));
	for (i = 0; i < m; i++)
		node_array[i].rem_cpus = 8;
	for (j = 0; j < n; j++) {
		job_array[j].min_nodes = 1;
		job_array[j].max_nodes = 4;
		job_array[j].priority = n*2 - j;
		job_array[j].min_cpus = (int) (j % 4 + 1) * 8;
	}
	i = solvewindow(m, n, 500, node_array, job_array);
}

