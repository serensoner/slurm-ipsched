#include <ilcplex/ilocplex.h>
#include "lpsched.h"
#include "src/plugins/select/lpconsres/select_lpconsres.h"

static void
populateproblem (IloModel model, IloNumVarArray x, IloRangeArray c,
			int resourceno, int nodeSize, int windowSize, 
			sched_nodeinfo_t *node_array, solver_job_list_t *job_list)
{
	IloEnv env = model.getEnv();
	IloObjective obj = IloMaximize(env);

	// Adding variables
	// x_ij : x[0 - nodeSize*windowSize]
	int varx = 0;
	for (int i = 0; i < nodeSize; i++)
		for (int j = 0; j < windowSize; j++) {
                        //stringstream ss;
                        //ss << "x" << i << j;
                        x.add(IloNumVar(env, 0.0, IloInfinity, ILOINT));	
                        x[(i * windowSize + j)].setName(ss.str().c_str());
                }	
		
	// s_j : x[nodeSize*windowSize :  (nodeSize + 1)*windowSize]
	for (int j = 0; j < windowSize; j++) {
		//stringstream ss;
		//ss << "s" << j;
		x.add(IloNumVar(env, 0.0, 1.0, ILOBOOL));
		x[nodeSize * windowSize + j].setName(ss.str().c_str());
	}
	
	// c_j : x[(nodeSize+1)*windowSize + (nodeSize+2)*windowSize]
	for (int j = 0; j < windowSize; j++) {
		//stringstream ss;
		//ss << "c" << j;
		x.add(IloNumVar(env, 0.0, 1.0, ILOFLOAT));
		x[(nodeSize + 1) * windowSize + j].setName(ss.str().c_str());
	}

	// t_ji : x[(nodeSize + 2) * windowSize : ((2 * nodeSize + 2)*windowSize]
	int varcounter = (nodeSize + 2) * windowSize;
	for (int j = 0; j < windowSize; j++) {
		for (int i = 0; i < nodeSize; i++) {
			//stringstream ss;
			//ss << "t" << j << i;
                        x.add(IloNumVar(env, 0.0, 1.0, ILOBOOL));
			x[varcounter + j * nodeSize + i].setName(ss.str().c_str());
		}
	}

	// sum over nodes (allocated cpu) should be equal to job's requested cpu
	// sum_j(x_ij) == r_j * s_j
	for (int j = 0; j < windowSize; j++) {
                IloExpr expr(env);
                for (int i = 0; i < nodeSize; i++)
                        expr += x[i * windowSize + j];
                model.add(expr == (int)job_list[j].min_cpus * x[nodeSize * windowSize + j]);
        }

        // sum over jobs for cpu should be available 
	// sum_j(x_ij) <= R_i
	for (int i = 0; i < nodeSize; i++) {
                c.add(IloRange(env, 0, (int)node_array[i].rem_cpus);
                for (int j = 0; j < windowSize; j++) 
                        c[i].setLinearCoef(x[(i * windowSize + j)], 1);
                //stringstream ss;
                //ss << "R" << i;
                c[i].setName(ss.str().c_str());
        }
	
        // sum over jobs for gpu should be available on each node
	// sum_j(t_ji * g_j) <= G_i
	for (int i = 0; i < nodeSize; i++) {
                c.add(IloRange(env, 0, (int)node_array[i].rem_gpus);
                for (int j = 0; j < windowSize; j++) 
                        c[i].setLinearCoef(x[varcounter + j * nodeSize + i], (int)job_array[j].gpu);
                //stringstream ss;
                //ss << "R" << i;
                c[i].setName(ss.str().c_str());
        }
	
	// gpu constraint here
        // t_ji is 1 if job j is allocated resourced in node i
	// t_ji = 1 if x_ij > 0
	for (int j = 0; j < windowSize; j++)  {
		for (int i = 0; i < nodeSize; i++) {
			IloExpr exp(env);
                        exp = (x[(i * windowSize + j)] > 0);
			model.add(exp == x[varcounter + j * nodeSize + i]);
		}
	}
     
        // this is not topologically aware, so cost is only selected nodes 
	// c_j constraint; sum_i(t_ji > 0) / (nodeSize + 1)
	for (int j = 0; j < windowSize; j++)  {
		IloExpr expr(env);
		for (int i = 0; i < nodeSize; i++)
                        expr += (x[varcounter + j * nodeSize + i]); 
		model.add( expr/ (nodeSize+1) == x[(nodeSize + 1) * windowSize + j]);
	}

	// min_nodes_j <= c_j * (nodeSize + 1) <= max_nodes_j
	for (int j = 0; j < windowSize; j++) {
		IloExpr exp(env);
		exp = (x[(nodeSize + 1) * windowSize + j] * nodeSize + 1);
		model.add(exp <= (int)job_array[j].max_nodes);
		model.add(exp >= (int)job_array[j].min_nodes);
	}

	for (int j = 0; j < windowSize; j++) {
		// p_j * s_j
		obj.setLinearCoef(x[nodeSize * windowSize + j], job_list[j].priority);
		// -p_j * c_j
		obj.setLinearCoef(x[(nodeSize + 1) * windowSize + j], -job_list[j].priority);
	}
	model.add(obj);
	model.add(c);
}

// resourceno is taken as 1, we are scheduling only for CPU!
void solve_allocation(int nodeSize, int windowSize, int timeout,
			sched_nodeinfo_t *node_array, solver_job_list_t *job_list)
{
	IloEnv   env;
	try {
		IloModel model(env);
		IloNumVarArray var(env);
		IloRangeArray con(env);

		int resourceno = 1;
		int currentTime = S->getTime();
		vector <float> p;
		vector <vector <int> > nodeLimit, r;
		
		populateproblem(model, var, con, resourceno, nodeSize, windowSize, node_array, job_list);

		IloCplex cplex(model);

                cplex.setParam(IloCplex::SimDisplay,  0);
                cplex.setParam(IloCplex::MIPDisplay,  0);
                cplex.setParam(IloCplex::SiftDisplay, 0);
                cplex.setParam(IloCplex::NetDisplay,  0);
                                         
                //stringstream filename;
                //filename << "model." << currentTime << ".lp";              
		//cplex.exportModel(filename.str().c_str());
                if ( !cplex.solve() ) {
			env.error() << "Failed to optimize LP" << endl;
			throw(-1);
		}

		IloNumArray vals(env);
		cplex.getValues(vals, var);
		//env.out() << "Solution status = " << cplex.getStatus() << endl;
		//env.out() << "Solution value  = " << cplex.getObjValue() << endl;
		//env.out() << "Values        = " << vals << endl;
		
		for (int j = 0; j < windowSize; j++) {
			if (cplex.getValue(var[nodeSize * windowSize + j]) > 0.99) { /* if the job is selected */
				solver_job_ptr = &job_list[j];
				solver_job_ptr->alloc_total = solver_job_ptr->min_cpus;
				solver_job_ptr->node_bitmap = (bitstr_t *) bit_alloc (node_record_count);
				solver_job_ptr->job_ptr->details->req_node_bitmap = (bitstr_t *) bit_alloc (node_record_count);
				solver_job_ptr->onnodes = (uint32_t *) xmalloc (sizeof(uint32_t)*node_record_count);
				for (int i = 0; i < node_record_count; i++)
					solver_job_ptr->onnodes[i] = cplex.getValue(var[(i * windowSize + j)]);
			}
		}
      	}
	catch (IloException& e) {
		cerr << "Concert exception caught: " << e << endl;
	}
	catch (...) {
		cerr << "Unknown exception caught at time " << S->getTime() << endl;
	}
	
	env.end();
} 

