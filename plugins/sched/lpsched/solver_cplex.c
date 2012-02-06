#include <solver_cplex.h>

static int populatebynonzero (CPXENVptr env, CPXLPptr lp, int m, int n, 
		int timeout, sched_nodeinfo_t *node_array, 
		solver_job_list_t *job_array);

static void free_and_null (char **ptr);

extern char* get_cplex_license_address(void);

extern int solve_allocation(int nodeSize, int windowSize, int timeout, 
			sched_nodeinfo_t *node_array, 
			solver_job_list_t *job_array)
{
	solver_job_list_t *sjob_ptr;
	struct job_details *job_det_ptr;
	int solstat;
	int n = windowSize, m = nodeSize;
	double objval;
	double *x = NULL;
	double *pi = NULL;
	double *slack = NULL;
	double *dj = NULL;

	CPXENVptr env = NULL;
	CPXLPptr lp = NULL;
	int status = 0;
	int i, j, k;
	int cur_numrows, cur_numcols;

	char envstr[256];
	sprintf(envstr,"ILOG_LICENSE_FILE=%s",get_cplex_license_address());
	if ( envstr != NULL ) {
		CPXputenv (envstr);
	}

	env = CPXopenCPLEX (&status);
	if ( env == NULL ) {
		char  errmsg[1024];
		fatal ("Could not open CPLEX environment.\n");
		CPXgeterrorstring (env, status, errmsg);
		fprintf (stderr, "%s", errmsg);
		goto TERMINATE;
	}

	status = CPXsetintparam (env, CPX_PARAM_SCRIND, CPX_ON);
	if (status) {
		fatal("Failure to turn on screen indicator, error %d.",status);
		goto TERMINATE;
	}

	status = CPXsetintparam (env, CPX_PARAM_DATACHECK, CPX_ON);
	if (status) {
		fatal("Failure to turn on data checking, error %d.", status);
		goto TERMINATE;
	}

	lp = CPXcreateprob (env, &status, "lpex1");

	if (lp == NULL) {
		fatal("Failed to create LP.");
		goto TERMINATE;
	}

	status = CPXsetdblparam(env,CPX_PARAM_TILIM,timeout);
	status = populatebynonzero (env, lp, m, n, timeout, 
					node_array, job_array);

	if ( status ) {
		fprintf (stderr, "Failed to populate problem.");
		goto TERMINATE;
	}

	status = CPXlpopt (env, lp);
	if ( status ) {
		fprintf (stderr, "Failed to optimize LP.");
		goto TERMINATE;
	}

	cur_numrows = CPXgetnumrows (env, lp);
	cur_numcols = CPXgetnumcols (env, lp);
	x = (double *) malloc (cur_numcols * sizeof(double));
	slack = (double *) malloc (cur_numrows * sizeof(double));
	dj = (double *) malloc (cur_numcols * sizeof(double));
	pi = (double *) malloc (cur_numrows * sizeof(double));

	if ( x == NULL ||
		slack == NULL ||
		dj    == NULL ||
		pi    == NULL   ) {
		status = CPXERR_NO_MEMORY;
		fatal ("Could not allocate memory for solution.");
		goto TERMINATE;
	}

	status = CPXsolution (env, lp, &solstat, &objval, x, pi, slack, dj);
	if ( status ) {
		fatal ("Failed to obtain solution.");
		goto TERMINATE;
	}

	for (j = 0; j < windowSize; j++) {
		if (x[j] > 0) {
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
			for (i = 0; i < nodeSize; i++) {
				k = (1 + i) * windowSize + j;
				if (x[k] > 0) {
					bit_set (sjob_ptr->node_bitmap, 
						(bitoff_t) (i));
					bit_set (job_det_ptr->req_node_bitmap, 
						(bitoff_t) (i));		
					node_array[i].rem_cpus -= x[k];
					node_array[i].rem_gpus -= sjob_ptr->gpu;
					sjob_ptr->onnodes[i] = x[k]; 
					job_det_ptr->req_node_layout[i] = 
						sjob_ptr->onnodes[i]; 
					sjob_ptr->alloc_total += x[k];
				}
			}
		} else
			job_array[j].alloc_total = 0;
	} 

TERMINATE:

	free_and_null ((char **) &x);
	free_and_null ((char **) &slack);
	free_and_null ((char **) &dj);
	free_and_null ((char **) &pi);

	if (lp != NULL) {
		status = CPXfreeprob (env, &lp);
		if (status) {
			fatal("CPXfreeprob failed, error code %d.", status);
		}
	}

	if (env != NULL) {
		status = CPXcloseCPLEX (&env);
		if (status) {
			char errmsg[1024];
			fatal("Could not close CPLEX environment.");
			CPXgeterrorstring (env, status, errmsg);
			fatal("%s", errmsg);
		}
	}     
	
	return (status);
}

static void free_and_null (char **ptr)
{
	if ( *ptr != NULL ) {
		free (*ptr);
		*ptr = NULL;
	}
}

inline int
populatebynonzero (CPXENVptr env, CPXLPptr lp, int m, int n, int timeout,
		sched_nodeinfo_t *node_array, solver_job_list_t *job_array)
{
	int NUMCOLS = n * (2 * m + 2);
	int NUMROWS = n * 3 + m * 2 + m * n * 2;
	int NUMNZ = 2 * n * (4 * m + 1); 
	int NZc = 0; /* nonzero counter */
	
	int status = 0;
	double *obj = NULL;
	obj = (double*)malloc(NUMCOLS * sizeof(double));
	double *lb = (double*)malloc(NUMCOLS * sizeof(double));
	double *ub = (double*)malloc(NUMCOLS * sizeof(double));
	double *rhs = (double*)malloc(NUMROWS * sizeof(double));
	char *sense = (char*)malloc(NUMROWS * sizeof(char));
	/*
	char **colname = (char**)malloc(NUMCOLS * sizeof(char*));
	char **rowname = (char**)malloc(NUMROWS * sizeof(char[10]));
	char str[10];
	*/
	int *rowlist = (int*)malloc(NUMNZ * sizeof(int));
	int *collist = (int*)malloc(NUMNZ * sizeof(int));
	double *vallist = (double*)malloc(NUMNZ * sizeof(double));
	int i, j, d;
	/*int rc = 0, cc = 0;*/

	CPXchgobjsen (env, lp, CPX_MAX);  /* Problem is maximization */

	/* row definitions */
	
	for (j = 0; j < n; j++) {
		sense[j] = 'E';
		rhs[j] = 0.0;
		/*debug3("cpueq row counter: %d, no: %d",rc++, j);*/
		/*sprintf(str,"CPUeq_%d",j+1);
		rowname[j] = str;*/
	}
	for (i = 0; i < m; i++) {
		sense[n + i] = 'L';
		rhs[n + i] = (double)(node_array[i].rem_cpus);
		/*sprintf(str,"NODEeq_%d",i+1);
		rowname[n + i] = str;*/
		/*debug3("nodeeq row counter: %d, no: %d",rc++, n+i);*/

		sense[n + m + i] = 'L';
		rhs[n + m + i] = (double)(node_array[i].rem_gpus);
		/*debug3("gpueq row counter: %d, no: %d",rc++, n+m+i);*/
		/*sprintf(str,"GPUeq_%d",i+1);
		rowname[n + m + i] = str;*/
	}

	for (j = 0; j < n; j++) {
		for (i = 0; i < m; i++) {
			d = n + 2 * m + j * m * 2 + i * 2;
			sense[d] = 'L';
			rhs[d] = 0.0;
			/*debug3("t_i_j_le row counter: %d, no: %d",
				rc++,d);*/
			/*sprintf(rowname[d],"t_%d_%d_LE",j+1,i+1);
			rowname[d] = str;*/
			sense[d + 1] = 'G';
			rhs[d + 1] = 0.0;
			/*debug3("t_i_j_ge row counter: %d, no: %d",
				rc++,d+1);*/
			/*sprintf(rowname[d+1],"t_%d_%d_GE",j+1,i+1);
			rowname[d+1] = str;*/
		}
	}

	for (j = 0; j < n; j++) {
		sense[n + 2 * m + n * m * 2 + j * 2] = 'L';
		rhs[n + 2 * m + n * m * 2 + j * 2] = 0.0;
		/*debug3("minmax row counter: %d, no: %d",
			rc++,n + 2 * m + n * m * 2 + j * 2);*/
		sense[n + 2 * m + n * m * 2 + j * 2 + 1] = 'G';
		rhs[n + 2 * m + n * m * 2 + j * 2 + 1] = 0.0;
		/*debug3("minmax row counter: %d, no: %d",
			rc++,n + 2 * m + n * m * 2 + j * 2 + 1);*/
		/*sprintf(str,"minmax_%d",j+1);
		rowname[n + 4 * m + n * m * 2 + j] = str;*/
	}

	status = CPXnewrows (env, lp, NUMROWS, rhs, sense, NULL, NULL);
	if ( status ) goto TERMINATE;

	/*debug3("ROWS: %d, column def starting",NUMROWS);*/
	/* column definitions */	
	for (j = 0; j < n; j++) {
		/*sprintf(str, "s_%d",j+1);
		colname[j] = str;*/
		lb[j] = 0.0;
		ub[j] = 1.0;
		/*debug3("s_j col counter: %d, no: %d",cc++,j);*/
	}

	/*debug3("defined s_j");*/
	for (i = 0; i < m; i++) {
		for (j = 0; j < n; j++) {
			/*sprintf(str, "x_%d_%d",i+1,j+1);
			colname[(i + 1) * n + j] = str;*/
			lb[(i + 1) * n + j] = 0.0;
			ub[(i + 1) * n + j] = CPX_INFBOUND;
			/*debug3("x_i_j col counter: %d, no: %d",
				cc++,(i+1)*n+j);*/
		}
	}
	/*debug3("defined x_i_j");*/

	for (j = 0; j < n; j++) {
		/*sprintf(str, "c_%d",j+1);
		colname[n * (m + 1) + j] = str;*/
		/* min_nodes_j <= c_j <= max_nodes_j */
		lb[n * (m + 1) + j] = (m + 1) * (job_array[j].min_nodes);
		ub[n * (m + 1) + j] = (m + 1) * (job_array[j].max_nodes);
		/*debug3("c_j col counter: %d, no: %d",cc++,n*(m+1)+j);*/
	}
	/*debug3("defined c_j");*/

	for (j = 0; j < n; j++) {
		for (i = 0; i < m; i++) {
			/*sprintf(str, "t_%d_%d",j+1,i+1);
			colname[n * (m + 2) + j * m + i] = str;*/
			lb[n * (m + 2) + j * m + i] = 0.0;
			ub[n * (m + 2) + j * m + i] = 1.0;
			/*debug3("t_i_j col counter: %d, no: %d",
				cc++,n*(m+2)+j*m+i);*/
		}
	}
	/*debug3("defined t_i_j");*/
	
	for (j = 0; j < NUMCOLS; j++) {
		obj[j] = 0;
	}
	for (j = 0; j < n; j++) {
		obj[j] = (job_array[j].priority);
		obj[n * (m + 1) + j] = (-1.0)*(job_array[j].priority);
	}

	status = CPXnewcols (env, lp, NUMCOLS, obj, lb, ub, NULL, NULL);
	if ( status )  goto TERMINATE;

	/*debug3("constraint coefficients");*/
	/* constraints */
	/* sum over nodes should be equal to job's required cpu */
	/* sum_i(x_ij) = r_j * s_j */
	for (j = 0; j < n; j++) {
		for (i = 0; i < m; i++) {
			rowlist[NZc] = j;			
			vallist[NZc] = 1.0;
			collist[NZc++] = (i + 1) * n + j;
			/*debug3("con1 collist %d",(i + 1) * n + j);*/
		}
		rowlist[NZc] = j;
		vallist[NZc] = (int)(-job_array[j].min_cpus);
		collist[NZc++] = j;
		/*debug3("con1 collist %d",j);*/
	}	
	/* nzc = m*n+n */

	/* sum over jobs for cpu should be available on nodes */
	/* sum_j(x_ij) <= R_i */
	for (i = 0; i < m; i++) {
		for (j = 0; j < n; j++) {
			/*debug3("con2 collist %d",(i + 1) * n + j);*/
			rowlist[NZc] = n + i;
			vallist[NZc] = 1.0;
			collist[NZc++] = (i + 1) * n + j;
		}
	}
	/* nzc = 2*m*n+n */

	/* sum over jobs for gpu should be available on nodes */
	/* sum_j(t_ji * g_j) <= G_i */
	for (i = 0; i < m; i++) {
		for (j = 0; j < n; j++) {
			rowlist[NZc] = n + m + i;
			vallist[NZc] = job_array[j].gpu;
			collist[NZc++] = n * (m + 2) + j * m + i;
			/*debug3("con3 collist %d",n * (m + 2) + j * m + i);*/
		}
	}
	/* nzc = 3*m*n+n */

	status = CPXchgcoeflist (env, lp, NZc, rowlist, collist, vallist);   
	if ( status )  goto TERMINATE;

TERMINATE:
	free_and_null ((char **) &obj);
	free_and_null ((char **) &lb);
	free_and_null ((char **) &ub);
	free_and_null ((char **) &rhs);
	free_and_null ((char **) &sense);
	return (status);
} 

