#include <ilcplex/cplex.h>

static void
	free_and_null (char **ptr);

typedef struct sched_nodeinfo {
	int rem_cpus;
	int rem_gpus;
} sched_nodeinfo_t;

typedef struct solver_job_list {
	int job_id;
	int min_nodes;
	int max_nodes;
	int min_cpus;
	int priority;
	int max_cpus;
	int alloc_total;
	int *onnodes;
	int gpu;
} solver_job_list_t;

extern int solve_allocation(int nodeSize, int windowSize, int timeout, 
			sched_nodeinfo_t *node_array, 
			solver_job_list_t *job_array);

int main() {
	int nodeSize = 1000;
	int windowSize = 500;
	int timeout = 0;
	int i;
	sched_nodeinfo_t *node_array = (sched_nodeinfo_t*)malloc(nodeSize * sizeof(sched_nodeinfo_t));
	solver_job_list_t *job_array = (solver_job_list_t*)malloc(windowSize * sizeof(solver_job_list_t));
	solver_job_list_t *solver_job_ptr;
	int minprio = 0;
	for (i = 0; i < nodeSize; i++) {
		node_array[i].rem_gpus = 0;
		node_array[i].rem_cpus = 16;
	}

	for (i = 0; i < windowSize; i++) {
		solver_job_ptr = &job_array[i];
		solver_job_ptr->min_nodes = 0;
		solver_job_ptr->max_nodes = 50000;
		solver_job_ptr->gpu = 0;
		solver_job_ptr->min_cpus = 16;
		solver_job_ptr->max_cpus = 16;
		solver_job_ptr->priority = 50000-i;
		if ((double)solver_job_ptr->priority < (double)minprio)
			minprio = (double)solver_job_ptr->priority;
		if ((int)solver_job_ptr->max_cpus < (int)solver_job_ptr->min_cpus)
			solver_job_ptr->max_cpus = solver_job_ptr->min_cpus;

	}
	solve_allocation(nodeSize, windowSize, timeout, node_array, job_array);
}


