/*
 * Iterative solver for heat distribution
 */

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "heat.h"

void usage(char *s) {
	fprintf(stderr, "Usage: %s <input file> [result file]\n\n", s);
}

inline void swap(double** a, double** b) {
	double *tmp = *a;
	*a = *b;
	*b = tmp;
}

int main(int argc, char *argv[]) {
	unsigned iter;
	FILE *infile, *resfile;
	char *resfilename;
	int myid, numprocs, my_x, my_y, block_size;
	int rows_per_proc, my_startrow;
	MPI_Status status, status2[2];

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);

	int last_proc = (myid == numprocs - 1);

	if (myid == 0) {

		// algorithmic parameters
		algoparam_t param;
		int np;
		double runtime, flop;
		double residual=0.0;

		// check arguments
		if( argc < 2 ) {
			usage( argv[0] );
			return 1;
		}

		// check input file
		if(!(infile=fopen(argv[1], "r"))) {
			fprintf(stderr, "\nError: Cannot open \"%s\" for reading.\n\n", argv[1]);
			usage(argv[0]);
			return 1;
		}

		// check result file
		resfilename = (argc>=3) ? argv[2] : "heat.ppm";

		if(!(resfile=fopen(resfilename, "w"))) {
			fprintf(stderr, "\nError: Cannot open \"%s\" for writing.\n\n", resfilename);
			usage(argv[0]);
			return 1;
		}

		// check input
		if(!read_input(infile, &param)) {
			fprintf(stderr, "\nError: Error parsing input file.\n\n");
			usage(argv[0]);
			return 1;
		}

		// check number of processes
		if((param.resolution % numprocs) != 0) {
			printf("[CRITICAL] number of processes dont split space evenly!\n");
			MPI_Abort(MPI_COMM_WORLD, -1);
		}

		// print parameters
		printf("MPI Num Procs     : %d\n", numprocs);
		print_params(&param);

		// set the visualization resolution
		param.u	 = 0;
		param.uhelp = 0;
		param.uvis  = 0;
		param.visres = param.resolution;

		if(!initialize(&param)) {
			fprintf(stderr, "Error in Solver initialization.\n\n");
			usage(argv[0]);
			return 1;
		}

		// full size (param.resolution are only the inner points)
		np = param.resolution + 2;

		// process parameters
		rows_per_proc = param.resolution / numprocs;
		my_startrow = rows_per_proc * myid;
	
		// starting time
		runtime = wtime();

		// send to workers the necessary data to perform computation
		MPI_Bcast(&param.maxiter,    1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&param.resolution, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&param.algorithm,  1, MPI_INT, 0, MPI_COMM_WORLD);
		for (int i=1; i<numprocs; i++) {
			if (i>0) {
				MPI_Send(&param.u[rows_per_proc*np*i],     (rows_per_proc+2)*np, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
				MPI_Send(&param.uhelp[rows_per_proc*np*i], (rows_per_proc+2)*np, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
			}
		}

		// print node information
		printf("[%d] working on rows %d to %d\n", myid, my_startrow, rows_per_proc*(myid+1)-1);

		iter = 0;

		while(1) {
			switch(param.algorithm) {
				case 0: // JACOBI
					residual = relax_jacobi(param.u, param.uhelp, rows_per_proc+2, np);
					swap(&param.u, &param.uhelp);
					// SWAP BOUNDARIES INFORMATION
					if(!last_proc) {
						// send data downstream
						MPI_Send(&param.u[(rows_per_proc)*np], np, MPI_DOUBLE, myid+1, iter, MPI_COMM_WORLD);
						// receive data upstream
						MPI_Recv(&param.u[(rows_per_proc+1)*np], np, MPI_DOUBLE, myid+1, iter, MPI_COMM_WORLD, &status);
					}
					break;
				case 1: // RED-BLACK - NOT WORKING WITH MPI!
					residual = relax_redblack(param.u, np, rows_per_proc+2);
					break;
				case 2: // GAUSS
					// SWAP BOUNDARIES INFORMATION
					residual = relax_gauss(param.u, rows_per_proc+2, np);
					if(!last_proc) {
						// send data downstream from current round
						MPI_Send(&param.u[(rows_per_proc)*np], np, MPI_DOUBLE, myid+1, iter, MPI_COMM_WORLD);
						// receive data upstream from current round
						MPI_Recv(&param.u[(rows_per_proc+1)*np], np, MPI_DOUBLE, myid+1, iter, MPI_COMM_WORLD, &status);
					}
					break;
			}

			iter++;

			// compute global residual
			MPI_Allreduce(MPI_IN_PLACE, &residual, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

			// solution good enough ?
			if (residual < 0.00005) break;

			// max. iteration reached ? (no limit with maxiter=0)
			if (param.maxiter>0 && iter>=param.maxiter) break;
		}

		// receive results from slaves
		for(int i=1; i<numprocs; i++) {
			MPI_Recv(&param.u[(rows_per_proc*i+1)*np], (rows_per_proc)*np, MPI_DOUBLE, i, param.maxiter+1, MPI_COMM_WORLD, &status);
		}

		// Flop count after iter iterations
		flop = iter * 11.0 * param.resolution * param.resolution;
		// stopping time
		runtime = wtime() - runtime;

		// for plot...
		coarsen(param.u, np, np, param.uvis, param.visres+2, param.visres+2);
		write_image(resfile, param.uvis, param.visres+2, param.visres+2);

		// print statistics
		fprintf(stdout, "[%d] finished computing with residual value = %f\n", myid, residual);
		fprintf(stdout, "[%d] Time: %04.3f ", myid, runtime);
		fprintf(stdout, "[%d] (%3.3f GFlop => %6.2f MFlop/s)\n", myid, flop/1000000000.0, flop/runtime/1000000);
		fprintf(stdout, "[%d] Convergence to residual=%f: %d iterations\n", myid, residual, iter);

		// cleanup
		finalize(&param);
		MPI_Finalize();

		return 0;

	} else {

		// receive information from master to perform computation locally
		int columns, rows, np;
		int rows_per_proc, my_startrow;
		int iter, maxiter;
		int algorithm;
		double residual;
		int last_proc = (myid == numprocs - 1);

		MPI_Bcast(&maxiter,    1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&columns, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&algorithm,  1, MPI_INT, 0, MPI_COMM_WORLD);

		// print slave information
		rows = columns;
		np = columns + 2;
		rows_per_proc = rows / numprocs;
		my_startrow = rows_per_proc * myid;
		
		// allocate memory for worker
		double * u     = calloc(sizeof(double),(rows_per_proc+2)*np);
		double * uhelp = calloc(sizeof(double),(rows_per_proc+2)*np);
		if((!u) || (!uhelp))	{
			fprintf(stderr, "Error: Cannot allocate memory\n");
			return 0;
		}
	
		// fill initial values for matrix with values received from master
		MPI_Recv(&u[0],     (rows_per_proc+2)*np, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
		MPI_Recv(&uhelp[0], (rows_per_proc+2)*np, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);

		// print node information
		printf("[%d] working on rows %d to %d\n", myid, my_startrow, rows_per_proc*(myid+1)-1);

		iter = 0;

		while(1) {
			switch(algorithm) {
				case 0: // JACOBI
					residual = relax_jacobi(u, uhelp, rows_per_proc+2, np);
					swap(&u, &uhelp);
					// SWAP BOUNDARIES INFORMATION
					// send data downstream 
					if (last_proc) MPI_Recv(&u[0], np, MPI_DOUBLE, myid-1, iter, MPI_COMM_WORLD, &status);
					else MPI_Sendrecv(&u[rows_per_proc*np], np, MPI_DOUBLE, myid+1, iter,
					                  &u[0],                np, MPI_DOUBLE, myid-1, iter,
					                  MPI_COMM_WORLD, &status);
					// send data upstream
					if (last_proc) MPI_Send(&u[np], np, MPI_DOUBLE, myid-1, iter, MPI_COMM_WORLD);
					else MPI_Sendrecv(&u[np],                   np, MPI_DOUBLE, myid-1, iter,
					                  &u[(rows_per_proc+1)*np], np, MPI_DOUBLE, myid+1, iter,
					                  MPI_COMM_WORLD, &status);
					break;

				case 1: // RED-BLACK - NOT WORKING WITH MPI!
					residual = relax_redblack(u, rows_per_proc+2, np);
					break;

				case 2: // GAUSS
					// SWAP BOUNDARIES INFORMATION
					// receive data downstream from current round (previous rows)
					MPI_Recv(&u[0], np, MPI_DOUBLE, myid-1, iter, MPI_COMM_WORLD, &status);
					// receive results upstream from last round
					if((!last_proc) && (iter != 0))
						MPI_Recv(&u[(rows_per_proc+1)*np], np, MPI_DOUBLE, myid+1, iter-1, MPI_COMM_WORLD, &status);
					residual = relax_gauss(u, rows_per_proc+2, np);
					// send data downstream from current round
					if(!last_proc)
						MPI_Send(&u[(rows_per_proc)*np], np, MPI_DOUBLE, myid+1, iter, MPI_COMM_WORLD);
					// send results upstream from current round
					MPI_Send(&u[np], np, MPI_DOUBLE, myid-1, iter, MPI_COMM_WORLD);
					break;
				}

			iter++;

			// compute global residual
			MPI_Allreduce(MPI_IN_PLACE, &residual, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

			// solution good enough ?
			if (residual < 0.00005) break;

			// max. iteration reached ? (no limit with maxiter=0)
			if (maxiter>0 && iter>=maxiter) break;
		}

		// send results to master
		MPI_Send(&u[np], (rows_per_proc)*np, MPI_DOUBLE, 0, maxiter+1, MPI_COMM_WORLD);

		// cleanup
		if(u) free(u);
		if(uhelp) free(uhelp);

		MPI_Finalize();
		exit(0);
	}
}
