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

#include "lpsched.h"
#include <ilcplex/cplex.h>

