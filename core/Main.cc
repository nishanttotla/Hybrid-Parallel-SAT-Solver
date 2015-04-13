/*****************************************************************************************[Main.cc]
Copyright (c) 2003-2006, Niklas Een, Niklas Sorensson
Copyright (c) 2007-2010, Niklas Sorensson

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
associated documentation files (the "Software"), to deal in the Software without restriction,
including without limitation the rights to use, copy, modify, merge, publish, distribute,
sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or
substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT
OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
**************************************************************************************************/
#include <omp.h>
#include <mpi.h>
#include <time.h>
#include <math.h>

#include <errno.h>

#include <signal.h>
#include <zlib.h>

#include "utils/System.h"
#include "utils/ParseUtils.h"
#include "utils/Options.h"
#include "core/Dimacs.h"
#include "core/Solver.h"


#define COMM_FACTOR .25


using namespace Minisat;

//=================================================================================================

int pow2(int in){
	if(in % 2 != 0)
		return 0;
	while(in % 2 == 0)
		in = in/2;
	return ((in > 1) ? 0 : 1);
}

void* master_comm_thread(int numprocs){
	int msg, restart_worker = 0, finished_count = 0, length = 0, count = 0;
	int* literals_msg;
	int* finished_workers = (int*) malloc(sizeof(int)*(numprocs + 1));
	printf("Master Comm Thread has started with %d workers\n", numprocs);
	
	while(finished_count != numprocs){
		//Wait for finished workers to send messages.
		MPI_Recv(&msg, 1, MPI_INT, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		printf("[MASTER_COMM_THREAD]::Received a restart message from Worker: %d\n",msg);
		finished_count++;
		printf("[MASTER_COMM_THREAD]:: %d Workers are now waiting for work\n",finished_count);
		//printf("[MASTER_COMM_THREAD]:: log2(finished_count + 1) is: %f and threshold is: %f\n",log2(finished_count + 1), COMM_FACTOR*numprocs);
		finished_workers[msg] = 1;
		
		//Check if we have a power of 2 number of solvers ready for restart.
		//Also need to make sure that this number is not too small. We don't want needless communication.
		//COMM_FACTOR is a tuning parameter which determines when restarts happen.
		if(pow2(finished_count + 1) == 1 && (finished_count + 1) >= COMM_FACTOR*(numprocs)){
			for(int i = 1; i <= numprocs; ++i){
				if(finished_workers[i] != 1){
					restart_worker = i;
					finished_workers[i] = 1;
					break;
				}
			}
			//Notify the worker who will seed restart.
			//The restart worker should see this message just before restarting its local SAT solver portfolio.
			//printf("[MASTER_COMM_THREAD]:: Number of workers restarting: %d\n", finished_count + 1);
			int literals_length = log(finished_count + 1)/log(2);
			printf("[MASTER_COMM_THREAD]:: Chose Worker %d as restart seed and now asking for %d vars\n",restart_worker, literals_length);
			MPI_Send(&literals_length, 1, MPI_INT, restart_worker, 1, MPI_COMM_WORLD);
			printf("[MASTER_COMM_THREAD]:: Send success(debug)\n");
			MPI_Recv(&length, 1, MPI_INT, restart_worker, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			printf("[MASTER_COMM_THREAD]:: Recv success(debug)\n");
			literals_msg = (int*) malloc(sizeof(int)*length);
			MPI_Recv(literals_msg, length, MPI_INT, restart_worker, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			printf("[MASTER_COMM_THREAD]:: Recv2 success(debug)\n");
			//Making sure that the restart happened successfully before resetting state variables. 
			//MPI_Recv(&msg, 1, MPI_INT, restart_worker, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			for(int i = 0; i <= numprocs; ++i){
				if(finished_workers[i] == 1 && i != restart_worker){
					printf("[MASTER_COMM_THREAD]:: Sending new guiding path to Worker %d\n", i);
					count++;
					MPI_Send(&count, 1, MPI_INT, i, 1, MPI_COMM_WORLD);
					MPI_Send(&length, 1, MPI_INT, i, 1, MPI_COMM_WORLD);
					MPI_Send(literals_msg, length, MPI_INT, i, 1, MPI_COMM_WORLD);
				}
			}
			printf("[MASTER_COMM_THREAD]:: Finished sending guiding paths and is checking for new workers\n");
			count = 0;
			finished_count = 0;
			for (int i = 0; i <= numprocs; i++)
				finished_workers[i] = 0;
		}
	}
	return NULL;
}


void printStats(Solver& solver)
{
    double cpu_time = cpuTime();
    double mem_used = memUsedPeak();
    printf("  ----------------------------------------         \n");
    printf("   winner:  more statistics           \n");
    printf("  ----------------------------------------         \n");
    printf("  restarts              : %"PRIu64"\n", solver.starts);
    printf("  conflicts             : %-12"PRIu64"   (%.0f /sec)\n", solver.conflicts   , solver.conflicts   /cpu_time);
    printf("  decisions             : %-12"PRIu64"   (%4.2f %% random) (%.0f /sec)\n", solver.decisions, (float)solver.rnd_decisions*100 / (float)solver.decisions, solver.decisions   /cpu_time);
    printf("  propagations          : %-12"PRIu64"   (%.0f /sec)\n", solver.propagations, solver.propagations/cpu_time);
    printf("  conflict literals     : %-12"PRIu64"   (%4.2f %% deleted)\n", solver.tot_literals, (solver.max_literals - solver.tot_literals)*100 / (double)solver.max_literals);
    if (mem_used != 0) printf("  Memory used           : %.2f MB\n", mem_used);
    printf("  CPU time              : %g s\n", cpu_time);
}


static Solver* solver;
// Terminate by notifying the solver and back out gracefully. This is mainly to have a test-case
// for this feature of the Solver as it may take longer than an immediate call to '_exit()'.
static void SIGINT_interrupt(int signum) { solver->interrupt(); }

// Note that '_exit()' rather than 'exit()' has to be used. The reason is that 'exit()' calls
// destructors and may cause deadlocks if a malloc/free function happens to be running (these
// functions are guarded by locks for multithreaded use).
static void SIGINT_exit(int signum) {
    printf("\n"); printf("*** INTERRUPTED ***\n");
    if (solver->verbosity > 0){
        printStats(*solver);
        printf("\n"); printf("*** INTERRUPTED ***\n"); }
    _exit(1); }


//=================================================================================================
// Main:


int main(int argc, char** argv)
{
	
    int numprocs, rank, rc;

	rc = MPI_Init(&argc, &argv);
	if (rc != MPI_SUCCESS){
		printf("Error starting MPI tasks. Terminating.\n");
		MPI_Abort(MPI_COMM_WORLD, rc);
	}

	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	printf("Number of Processors = %d My rank = %d\n",numprocs, rank);
	lbool result;
	try {
        setUsageHelp("USAGE: %s [options] <input-file> <result-output-file>\n\n  where input may be either in plain or gzipped DIMACS.\n");
        // printf("This is MiniSat 2.0 beta\n");
        
#if defined(__linux__)
        fpu_control_t oldcw, newcw;
        _FPU_GETCW(oldcw); newcw = (oldcw & ~_FPU_EXTENDED) | _FPU_DOUBLE; _FPU_SETCW(newcw);
        printf("WARNING: for repeatability, setting FPU to use double precision\n");
#endif
        // Extra options:
        //
        IntOption    verb   ("MAIN", "verb",   "Verbosity level (0=silent, 1=some, 2=more).", 1, IntRange(0, 2));
        IntOption    cpu_lim("MAIN", "cpu-lim","Limit on CPU time allowed in seconds.\n", INT32_MAX, IntRange(0, INT32_MAX));
        IntOption    mem_lim("MAIN", "mem-lim","Limit on memory usage in megabytes.\n", INT32_MAX, IntRange(0, INT32_MAX));

	IntOption    ncores ("MAIN", "ncores","# threads.\n", 4, IntRange(1, sysconf( _SC_NPROCESSORS_ONLN )));
	IntOption    limitEx("MAIN", "limitEx","Limit size clause exchange.\n", 10, IntRange(0, INT32_MAX));
    IntOption    det    ("MAIN", "det","determenistic mode (0=non deterministic, 1=deterministic static, 2=deterministic dynamic.\n", 0, IntRange(0, 2));
	IntOption    ctrl   ("MAIN", "ctrl","Dynamic control clause sharing with 2 modes.\n", 0, IntRange(0, 2));
	
	time_t start_time = time(NULL);
	parseOptions(argc, argv, true);

	double initial_time = cpuTime();

		
	int nbThreads   = ncores;
	int limitExport = limitEx;	
	Cooperation* coop= new Cooperation(nbThreads, limitExport);
	
	coop->ctrl = ctrl;
	coop->deterministic_mode = det;

	//printf("Diversifying portfolio...\n");
	//coop->diversify(nbThreads);
	coop->printPortfolio(nbThreads);
	
	for(int t = 0; t < nbThreads; t++){
	  coop->solvers[t].threadId = t;
	  coop->solvers[t].verbosity = verb;
	  coop->solvers[t].deterministic_mode = det;
	}

	
	if(rank == 0){
		//Launch comm thread
		master_comm_thread(numprocs-1);
	}

	printf(" -----------------------------------------------------------------------------------------------------------------------\n");
	printf("|                                 manysat2.0    %i thread(s) on %i core(s)                                                |\n", coop->nbThreads, sysconf( _SC_NPROCESSORS_ONLN )); 
	printf(" -----------------------------------------------------------------------------------------------------------------------\n");



		
        // Use signal handlers that forcibly quit until the solver will be able to respond to
        // interrupts:
        signal(SIGINT, SIGINT_exit);
        signal(SIGXCPU,SIGINT_exit);

        // Set limit on CPU-time:
        if (cpu_lim != INT32_MAX){
            rlimit rl;
            getrlimit(RLIMIT_CPU, &rl);
            if (rl.rlim_max == RLIM_INFINITY || (rlim_t)cpu_lim < rl.rlim_max){
                rl.rlim_cur = cpu_lim;
                if (setrlimit(RLIMIT_CPU, &rl) == -1)
                    printf("WARNING! Could not set resource limit: CPU-time.\n");
            } }

        // Set limit on virtual memory:
        if (mem_lim != INT32_MAX){
            rlim_t new_mem_lim = (rlim_t)mem_lim * 1024*1024;
            rlimit rl;
            getrlimit(RLIMIT_AS, &rl);
            if (rl.rlim_max == RLIM_INFINITY || new_mem_lim < rl.rlim_max){
                rl.rlim_cur = new_mem_lim;
                if (setrlimit(RLIMIT_AS, &rl) == -1)
                    printf("WARNING! Could not set resource limit: Virtual memory.\n");
            } }
        
        if (argc == 1)
            printf("Reading from standard input... Use '--help' for help.\n");
        
        gzFile in = (argc == 1) ? gzdopen(0, "rb") : gzopen(argv[1], "rb");
        if (in == NULL)
            printf("ERROR! Could not open file: %s\n", argc == 1 ? "<stdin>" : argv[1]), exit(1);
	
	

        if (coop->solvers[0].verbosity > 0){
            printf(" ===============================================[ Problem Statistics ]==================================================\n");
            printf("|                                                                                                                       |\n"); }
        
        //parse_DIMACS(in, &S);
        parse_DIMACS(in, coop);
		
		
		vec<Lit> lits;
		int i = numprocs-2;
		int numvars = 0;
		while (i > 0){
			i >>= 1;
			numvars++;
		}
		printf("Number of variables being guided %d\n",numvars);
		int id = rank;
		
		// default for debugging
		int staticvars[15] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
		
		//TLB_Module (no restarts)
		//int staticvars[15] = {1, 2203040, 2203039, 2203038, 2203037, 2203036, 2203035, 2203034, 2203033, 2203032, 2203031, 2203030, 2203029, 2203028, 2203027};
		// int staticvars[15] = {1, 2203040, 2203039, 2203038, 2203037, 2203036, 2203035, 2203034, 2203033, 2203032, 2203031, 2203030, 2203029, 2203028, 2203027};
		// int staticvars[15] = {142580, 145481, 142830, 142499, 134533, 143205, 125587, 136340, 160506, 43326, 145106, 142955, 219758, 124466, 124736};
		//
		
		//gss (no restarts)
		//int staticvars[15] = {1, 31821 ,31820, 31819, 31818 ,31817, 31816, 31815, 31814, 31813, 31812, 31811, 31810, 31809, 31808};
		//
		
		//f10bidw (no restarts)
		// int staticvars[15] = {1, 824632,824629,824627,824625,824623,824621,824619,824617,824615,824613,824611,824609,824607,824605};
		//

		//f10nidw (same category of instance as f10bidw)
		// int staticvars[15] = {69963, 69887, 69287, 71097, 42922, 69735, 32762, 71095, 152714, 43841, 42432, 69773, 54533, 25675, 31915};
		
		//simon03 (no restarts)
		// int staticvars[15] = {2593, 5161, 5132, 8632, 1249, 11213, 7204, 4745, 12308, 10745, 4246, 808, 527, 11340, 10765};
		//
		
		//uti (no restarts)
		//int staticvars[15] = {1, 259597 ,25959, 259595,259594,259593,259592,259591,259590,259589,259588,259587,259586,259585,259584};
		//
		
		//partial-10-11-s
		//int staticvars[15] = {13, 378911, 378908, 378907, 378905, 378903, 378901, 378899, 378897, 378895, 378893, 378891, 378889, 378887, 378885};
		//

		//partial-10-19-s
		//int staticvars[15] = {122919, 274001, 191079, 13, 661591, 661589, 661587, 661585, 661583, 661581, 661579, 661577, 661575, 661573, 661571};
		//

		//dated-5-17-u
		//int staticvars[15] = {13, 366927, 366925, 366923, 366921, 366919, 366917, 366915, 366913, 366911, 366909, 366907, 366905, 366903, 366901};
		//

		//partial-10-19-u
		//int staticvars[15] = {13, 660743, 660741, 660739, 660737, 660735, 660733, 660731, 660729, 660727, 660725, 660723, 660721, 660719, 660717};
		//
		
		while(numvars != 0 && numprocs != 1){
			int var = staticvars[ numvars - 1] - 1;
			int initvar = var;
			printf("MPI_RANK[%d]::setting variable: %d = %d\n", rank, staticvars[numvars - 1], (id % 2));
			for(int t = 0; t < coop->nbThreads; t++){
				var = initvar;
				while (var >= coop->solvers[t].nVars()) coop->solvers[t].newVar();
			}
			lits.push( ((id % 2) == 0) ? mkLit(var) : ~mkLit(var) );
			id = id >> 1;
			numvars--;
			for(int t = 0; t < coop->nbThreads;t++){
				coop->solvers[t].addClause(lits);
			}
			lits.clear();
		}
		
		gzclose(in);
        FILE* res = (argc >= 3) ? fopen(argv[2], "wb") : NULL;
        
        if (coop->solvers[0].verbosity > 0){
	  printf("|  Number of cores:      %12d                                                                                   |\n", coop->nbThreads); 
	  printf("|  Number of variables:  %12d                                                                                   |\n", coop->solvers[0].nVars());
	  printf("|  Number of clauses:    %12d                                                                                   |\n", coop->solvers[0].nClauses()); 

}
        
        double parsed_time = cpuTime();
        if (coop->solvers[0].verbosity > 0){
            printf("|  Parse time:           %12.2f s                                                                                 |\n", parsed_time - initial_time);
            printf("|                                                                                                                       |\n"); }
 

        // Change to signal-handlers that will only notify the solver and allow it to terminate
        // voluntarily:
        signal(SIGINT, SIGINT_interrupt);
        signal(SIGXCPU,SIGINT_interrupt);
       
 
        
    vec<Lit> dummy;	
	int winner = 0;
	lbool ret;
	int path_config;
	int* literals_msg;
	int restart = 0;
restart:

	if (coop->solvers[0].verbosity > 0){
	  printf(" ==============================================[ Search Statistics ]====================================================\n");
	  printf("|   Thread  | Conflicts |          ORIGINAL         |          LEARNT          |  exported  / imported      |  Progress |\n");
	  printf("|    (id)   |           |    Vars  Clauses Literals |    Limit  Clauses Lit/Cl |  (clauses by restart)      |           |\n");
	  printf(" =======================================================================================================================\n");
    }
	if (restart)
		printf("Worker[%d] just restarted and is ready to call the solver\n",rank);
	if (!coop->solvers[0].simplify()){
	  if (res != NULL) fprintf(res, "UNSAT\n"), fclose(res);
	  if (coop->solvers[0].verbosity > 0){
	    printf("========================================================================================================================\n");
	    printf("Solved by unit propagation\n");
	    printStats(coop->solvers[0]);
	    printf("\n"); }
	  	printf("Worker %d: UNSATISFIABLE\n", rank);
	  goto restart_comm;
	  //MPI_Finalize();
	  //exit(20);
    }
	// launch threads in Parallel 	
	omp_set_num_threads(nbThreads);
#pragma omp parallel
	{
	  int t = omp_get_thread_num();
	  coop->start = true; 
	  ret = coop->solvers[t].solveLimited(dummy, coop);
	}

	printf("Worker %d just finished solveLimited()\n", rank);
	
	// select winner threads with respect to deterministic mode 
	for(int t = 0; t < coop->nThreads(); t++)
	  if(coop->answer(t) != l_Undef){
	    winner = t;
	    result = coop->answer(t);
	    break;
	  }
	coop->printStats(coop->solvers[winner].threadId);
	printStats(coop->solvers[winner]);
	if(result == l_True){
		printf("SATISFIABLE\n");
		MPI_Abort(MPI_COMM_WORLD, 0);
	}
	printf(result == l_False ? "Worker %d UNSATISFIABLE\n" : "Worker %d INDETERMINATE\n", rank);
restart_comm:
	
	delete coop;
	coop = new Cooperation(nbThreads, limitExport);

	coop->ctrl = ctrl;
	coop->deterministic_mode = det;

	//printf("Diversifying portfolio...\n");
	//coop->diversify(nbThreads);
	coop->printPortfolio(nbThreads);
	
	for(int t = 0; t < nbThreads; t++){
	  coop->solvers[t].threadId = t;
	  coop->solvers[t].verbosity = verb;
	  coop->solvers[t].deterministic_mode = det;
	}

	in = (argc == 1) ? gzdopen(0, "rb") : gzopen(argv[1], "rb");
	if (in == NULL)
		printf("ERROR! Could not open file: %s\n", argc == 1 ? "<stdin>" : argv[1]), exit(1);

	//if (coop->solvers[0].verbosity > 0){
		printf(" ===============================================[ Worker %d Parsing Input File and Restarting ]==================================================\n", rank);
		printf("|                                                                                                                       |\n"); 
	//}
	
	//parse_DIMACS(in, &S);
	parse_DIMACS(in, coop);
	gzclose(in);
	//Notify Master_comm_thread that this process finished and ready for potential restart.
	MPI_Send(&rank, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
	int length;
	//need to clear learned clauses
	//need to clear all temporary variables
	//restart solver with new literals
	printf("Worker %d waiting for path_config and new guiding path\n", rank);
	MPI_Recv(&path_config, 1, MPI_INT, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	MPI_Recv(&length, 1, MPI_INT, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	printf("Worker %d received path_config: %d and length: %d\n",rank, path_config, length);
	if(length == -1)
		goto end;
	literals_msg = (int*) malloc(sizeof(int)*length);
	MPI_Recv(literals_msg, length, MPI_INT, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	
	printf("WORKER[%d]::Is ready to restart with new guiding path\n",rank);
	//Setup new literals.
	while(length != 0 && numprocs != 1){
		int var = literals_msg[ length - 1];
		int initvar = var;
		printf("MPI_RANK[%d]::setting variable: %d = %d\n", rank, literals_msg[length - 1], (path_config % 2));
		for(int t = 0; t < coop->nbThreads; t++){
			var = initvar;
			while (var >= coop->solvers[t].nVars()) coop->solvers[t].newVar();
		}
		lits.push( ((path_config % 2) == 0) ? mkLit(var) : ~mkLit(var) );
		path_config = path_config >> 1;
		length--;
		for(int t = 0; t < coop->nbThreads;t++){
			coop->solvers[t].addClause(lits);
		}
		lits.clear();
	}
	restart = 1;
	goto restart;

end:
	time_t end_time = time(NULL);
	double elapsed = difftime(end_time, start_time);
	printf("MPI_RANK[%d] took %f seconds to finish\n",rank, elapsed);

	if (coop->solvers[winner].verbosity > 0){
	  //printStats(coop->solvers[winner]);
            printf("\n"); }
        if(result == l_True){
			printf("SATISFIABLE\n");
			MPI_Abort(MPI_COMM_WORLD, 0);
			//goto stop;
		}
		
		printf(result == l_False ? "UNSATISFIABLE\n" : "INDETERMINATE\n");
        if (res != NULL){
            if (result == l_True){
	      fprintf(res, "SAT\n");
                for (int i = 0; i < coop->solvers[winner].nVars(); i++)
                    if (coop->solvers[winner].model[i] != l_Undef)
                        fprintf(res, "%s%s%d", (i==0)?"":" ", (coop->solvers[winner].model[i]==l_True)?"":"-", i+1);
                fprintf(res, " 0\n");
            	goto stop;
			}else if (result == l_False)
                fprintf(res, "UNSAT\n");
            else
                fprintf(res, "INDET\n");
            fclose(res);
        }
        
	if(rank == 0){
		void* tmp;
		//pthread_join(comm_thread, &tmp);
	}

    
	} catch (OutOfMemoryException&){
        printf("===============================================================================\n");
        printf("INDETERMINATE\n");
		goto stop;
    }
stop:
	printf("MPI_RANK[%d] is exiting\n",rank);
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
	
	
	#ifdef NDEBUG
			exit(result == l_True ? 10 : result == l_False ? 20 : 0);     // (faster than "return", which will invoke the destructor for 'Solver')
	#else
			return (result == l_True ? 10 : result == l_False ? 20 : 0);
	#endif
}
