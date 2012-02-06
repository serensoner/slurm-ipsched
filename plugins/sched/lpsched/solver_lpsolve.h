/*****************************************************************************\
 *  lpsched.c - linear programming based scheduler plugin.
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

#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#include "slurm/slurm.h"
#include "slurm/slurm_errno.h"

#include "src/common/list.h"
#include "src/common/macros.h"
#include "src/common/node_select.h"
#include "src/common/parse_time.h"
#include "src/common/slurm_protocol_api.h"
#include "src/common/xmalloc.h"
#include "src/common/xstring.h"

#include "src/slurmctld/acct_policy.h"
#include "src/slurmctld/front_end.h"
#include "src/slurmctld/job_scheduler.h"
#include "src/slurmctld/licenses.h"
#include "src/slurmctld/locks.h"
#include "src/slurmctld/node_scheduler.h"
#include "src/slurmctld/preempt.h"
#include "src/slurmctld/reservation.h"
#include "src/slurmctld/slurmctld.h"
#include "src/slurmctld/srun_comm.h"

#include "lp_lib.h"
#include "src/plugins/select/lpconsres/select_lpconsres.h"
#include "lpsched.h"

