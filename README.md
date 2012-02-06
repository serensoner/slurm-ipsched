lpsched - lp programming based co-allocation scheduler plugin.

Developed by Seren Soner and Can Ozturan from Bogazici University, Istanbul, Turkey as a part of PRACE WP9.2

Contact seren.soner@gmail.com, ozturaca@boun.edu.tr if necessary

- IMPORTANT NOTE: sched currently works with only CR_CPU with GresTypes=gpu (optional)
  nodes should be SHARED, overcommit is not allowed.

--------------------------------------------------------
The main idea is the following:

- from the queue, do not start the first job, but take a set of jobs 
- test all of them at once to see selection of which of these jobs 
  would yield the best result

--------------------------------------------------------
TODO

- (week of Dec 5-9) - fatal error in initalization for sched/lpsched if select/lpconsres is not selected and vice versa.
		    also fatal error is required when other SelectParameters are used.
- (week of Dec 5-9) - clean the code -- too much debugging comments both useful and useless.
- (week of Dec 5-9) - need a new makefile to compile the code with lpsolve dynamically. currently it uses static liblpsolve55.a for i386 systems.
- after the paper is published, a detailed readme, explaining the logic behind solve_allocation().

--------------------------------------------------------
CHANGES

- Dec 5, 2011: gpu support added. for this, common/gres.c was modified and some new functions were added. it returns the job's gpu request and node's availab
le gpu count to the lpsched scheduler. 
- Dec 4, 2011: submitted jobs are not scheduled right away currently. they should first go through solve_allocation(). this required a change in the core cod
e in the job_mgr.c line 2852. It does not call select_nodes() anymore.

