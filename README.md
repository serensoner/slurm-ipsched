ipsched - ip programming based co-allocation scheduler plugin for SLURM.

Developed by Seren Soner and Can Ozturan from Bogazici University, Istanbul, Turkey as a part of PRACE WP9.2

Contact seren.soner@gmail.com, ozturaca@boun.edu.tr if necessary

- IMPORTANT NOTE: sched currently works with only CR_CPU with GresTypes=gpu (optional)
  nodes should be SHARED, overcommit is not allowed.

- USE CPLEX WHENEVER AVAILABLE. In order to use cplex, first go to the plugins/sched/lpsched directory, copy Makefile.am.cplex to Makefile.am. If you wish to use lpsolve, copy Makefile.am.lpsolve to Makefile.am. Afterwards, run autogen.sh, ./configure, and make install. If CPLEX is selected, the license address of cplex should be given in slurm.conf in the following way: 

	`SchedulerParameters=cplex_lic=/share/apps/ILOG/CPLEX_Studio/AcademicResearch122/licenses/access.ilm`

- max_job_count can also be set to limit the number of jobs handled in every single ILP solve of the scheduler. By default, it is 200.

--------------------------------------------------------
The main idea is the following:

- from the queue, do not start the first job, but take a set of jobs 
- test all of them at once to see selection of which of these jobs 
  would yield the best result

--------------------------------------------------------
IMPORTANT NOTE ABOUT PATCHES:

- patches on gres.c and gres.h add two functions, which return number of gres types defined, and the number of gpu's a specific job requests. also, one function is modified to return number of gpu's available in that node, instead of not returning anything.

- patch on job_scheduler.c changes the job scheduler so that, when a job is submitted, its not immediately considered by the job scheduler; but waits for the next run of IP solve to be allocated.

- patch on makefile.am files allow the makefiles to be generated automatically.


--------------------------------------------------------
CHANGES

- Feb 6, 2011: Restructuring.
- Jan 31, 2011: Bug fixes.
- Jan 10, 2011: cplex and lpsolve can now be compiled dynamically.
- Dec 5, 2011: gpu support added. for this, common/gres.c was modified and some new functions were added. it returns the job's gpu request and node's availab
le gpu count to the lpsched scheduler. 
- Dec 4, 2011: submitted jobs are not scheduled right away currently. they should first go through solve_allocation(). this required a change in the core cod
e in the job_mgr.c line 2852. It does not call select_nodes() anymore.

