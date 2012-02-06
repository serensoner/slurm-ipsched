/*****************************************************************************\
 *  lpsched.c - simple backfill scheduler plugin.
 *
 *  Developed by Seren Soner and Can Ozturan 
 *  from Bogazici University, Istanbul, Turkey
 *  as a part of PRACE WP9.2
 *  Contact seren.soner@gmail.com, ozturaca@boun.edu.tr if necessary
 
 *  Used backfill plugin as base code.

 *  For the theory behind the plugin, see the manuscript
 *  doi: ....
\*****************************************************************************/

#ifndef _SLURM_LPSCHED_H
#define _SLURM_LPSCHED_H

typedef struct node_space_map {
	time_t begin_time;
	time_t end_time;
	bitstr_t *avail_bitmap;
	int next;	/* next record, by time, zero termination */
} node_space_map_t;

typedef struct select_nodeinfo {
	uint16_t magic;		/* magic number */
	uint16_t alloc_cpus;
	uint16_t rem_cpus;
	uint16_t alloc_mem;
	uint16_t rem_mem;
	uint16_t rem_gpus;
} select_node_info_t ;

typedef struct sched_nodeinfo {
	uint16_t rem_cpus;
	uint16_t rem_gpus;
} sched_nodeinfo_t;

extern char* get_cplex_license_address(void);

/* lpsched_agent - detached thread periodically attempts to backfill jobs */
extern void *lpsched_agent(void *args);

/* Terminate lpsched_agent */
extern void stop_lpsched_agent(void);

/* Note that slurm.conf has changed */
extern void lpsched_reconfig(void);

typedef struct solver_job_list {
	struct job_record *job_ptr;
	uint32_t job_id;
	uint32_t min_nodes;
	uint32_t max_nodes;
	uint32_t min_cpus;
	uint32_t priority;
	uint32_t max_cpus;
	uint32_t alloc_total;
	uint32_t *onnodes;
	bitstr_t *node_bitmap;
	uint32_t gpu;
} solver_job_list_t;

struct sched_nodeinfo * _print_nodes_inlpsched();
#endif	/* _SLURM_LPSCHED_H */
