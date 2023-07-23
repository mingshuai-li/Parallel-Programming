#include <stdio.h>
#include <stdlib.h>

#include "input.h"
#include "heat.h"
#include "timing.h"
#include <mpi.h>
#include <omp.h>

void usage(char *s) {
	fprintf(stderr, "Usage: %s <input file> <positive mpi topology rows> <positive mpi topology columns> [result file]\n\n", s);
}

int main(int argc, char *argv[]){
	
	// check arguments
	if (argc < 4) {
		usage(argv[0]);
		return 1;
	}	

	int i, j, k;
	FILE *infile, *resfile;
	char *resfilename;
	int iter;

    double* time;
	double localResidual = 0.0;
	double globalResidual = 0.0;

	int size, rank;
	MPI_Status status;
    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    MPI_Comm_size (MPI_COMM_WORLD, &size);
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);

	int mpiTopoRows = atoi(argv[2]); // the number of rows in the MPI virtual topology
	int mpiTopoColumns = atoi(argv[3]); // // the number of columns in the MPI virtual topology

	// check arguments
    if (mpiTopoRows<=0 || mpiTopoColumns<=0) {
		usage(argv[0]);
		return 1;
    }	
	
	MPI_Comm gridComm2d;
	int dims[2] = {mpiTopoRows, mpiTopoColumns};
	int coordinates[2] = {0, 0};
	int periods[2] = {0, 0};
	int reorder = 0;	
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, reorder, &gridComm2d);
    MPI_Cart_coords(gridComm2d, rank, 2, coordinates);

	// get ranks of neighbors
	enum directions {north, south, west, east};
	int neighborsRanks[4];
	MPI_Cart_shift(gridComm2d, 0, 1, &neighborsRanks[north], &neighborsRanks[south]);
	MPI_Cart_shift(gridComm2d, 1, 1, &neighborsRanks[west], &neighborsRanks[east]);
	
	// printf("I am rank %d, my coordinates are %d and %d, "
    //    "my north neighbor is rank %d, "
    //    "my south neighbor is rank %d, "
    //    "my west neighbor is rank %d, "
    //    "my east neighbor is rank %d. \n", 
    //    rank, coordinates[0], coordinates[1],
    //    neighborsRanks[north], neighborsRanks[south], 
    //    neighborsRanks[west], neighborsRanks[east]);

    algoparam_t param;
	param.visres = 100;

	// check input file
	if (!(infile = fopen(argv[1], "r"))) {
		fprintf(stderr, "\nError: Cannot open \"%s\" for reading.\n\n", argv[1]);
		usage(argv[0]);
		return 1;
	}

	// check result file
	resfilename = (argc >= 5) ? argv[4] : "heat.ppm";

	if (!(resfile = fopen(resfilename, "w"))) {
		fprintf(stderr, "\nError: Cannot open \"%s\" for writing.\n\n", resfilename);
		usage(argv[0]);
		return 1;
	}

	// check input
	if (!read_input(infile, &param)) {
		fprintf(stderr, "\nError: Error parsing input file.\n\n");
		usage(argv[0]);
		return 1;
	}

	// check algorithm (communication)
    if (param.algorithm < 0 || param.algorithm > 3) {
        fprintf(stderr, "\nError: Error specifying algorithm (communication).\n\n");
		usage(argv[0]);
		return 1;
    }

	if (rank == 0) {
		printf("\nMPI topology configuration: %d rows * %d columns. \n\n", mpiTopoRows, mpiTopoColumns);
		print_params(&param);
		time = (double *) calloc(sizeof(double), (int) (param.max_res - param.initial_res + param.res_step_size) / param.res_step_size);
	}

	int exp_number = 0;
	int remaining = 0;

	for (param.act_res = param.initial_res; param.act_res <= param.max_res; param.act_res = param.act_res + param.res_step_size) {
		
		// horizontally divide the square (without boundary), 
		// start/endRowIndex marks the location of the process in the global square with boundary
		// remaining is distributed evenly at the first <remaining> rows of processes
		param.numberOfRows = param.act_res / mpiTopoRows;
		remaining = param.act_res % mpiTopoRows;
		if(coordinates[0] < remaining){
			param.numberOfRows++;
			param.startRowIndex = coordinates[0] * param.numberOfRows + 1;
		}
		else{
			param.startRowIndex = coordinates[0] * param.numberOfRows + 1 + remaining;
		}
		param.endRowIndex = param.startRowIndex + param.numberOfRows - 1;

		// vertically divide the square (without boundary), 
		// start/endColumnIndex marks the location of the process in the global square with boundary
		// remaining is distributed evenly at the first <remaining> columns of processes
		param.numberOfColumns = param.act_res / mpiTopoColumns;
		remaining = param.act_res % mpiTopoColumns;
		if(coordinates[1] < remaining){
			param.numberOfColumns++;
			param.startColumnIndex = coordinates[1] * param.numberOfColumns + 1;
		}
		else{
			param.startColumnIndex = coordinates[1] * param.numberOfColumns + 1 + remaining;
		}
		param.endColumnIndex = param.startColumnIndex + param.numberOfColumns - 1;

		param.numberOfRowsIncludingGhost = param.numberOfRows + 2;
		param.numberOfColumnsIncludingGhost = param.numberOfColumns + 2;
		
	// 	printf("I am rank %d, my coordinates are %d and %d, "
    //    "I have %d rows, and %d rows with 2 ghost, "
    //    "Globally, I start from row %d to row %d, "
    //    "I have %d columns, and %d columns with 2 ghost, "
    //    "Globally, I start from column %d to column %d. \n", 
    //    rank, coordinates[0], coordinates[1],
    //    param.numberOfRows, param.numberOfRowsIncludingGhost, 
	//    param.startRowIndex, param.endRowIndex, 
	//    param.numberOfColumns, param.numberOfColumnsIncludingGhost, 
	//    param.startColumnIndex, param.endColumnIndex);
       	
		if (!initialize(&param)) {
			fprintf(stderr, "Error in Jacobi initialization.\n\n");
			usage(argv[0]);
		}

		// if(rank==0){
		// 	for(i = 0; i<param.numberOfRowsIncludingGhost;i++){
		// 		for(j = 0; j<param.numberOfColumnsIncludingGhost;j++){
		// 			printf("%f ", param.u[i*param.numberOfColumnsIncludingGhost+j]);
		// 		}
		// 		printf("\n");
		// 	}
		// }
		
		// west and east ghost columns are not contiguous and need tmp columns for communication
		double* sendWestColumn = (double*)malloc( sizeof(double)* param.numberOfRows);
		double* recvWestColumn = (double*)malloc( sizeof(double)* param.numberOfRows);
		double* sendEastColumn = (double*)malloc( sizeof(double)* param.numberOfRows);
		double* recvEastColumn = (double*)malloc( sizeof(double)* param.numberOfRows);
		

		MPI_Barrier(gridComm2d);
		if (rank==0) time[exp_number] = wtime();

		switch(param.algorithm){

		// blocking communication
		case 0:
		for (iter = 0; iter < param.maxiter; iter++) {
			localResidual = relax_jacobi_blocking(&(param.u), &(param.uhelp), param.numberOfColumnsIncludingGhost, param.numberOfRowsIncludingGhost, 
												  &gridComm2d, neighborsRanks, &sendWestColumn, &recvWestColumn, &sendEastColumn, &recvEastColumn);
			MPI_Reduce(&localResidual, &globalResidual, 1, MPI_DOUBLE, MPI_SUM, 0, gridComm2d);
		}
		break;

		// non-blocking communication
		case 1:
		for (iter = 0; iter < param.maxiter; iter++) {
			localResidual = relax_jacobi_nonblocking(&(param.u), &(param.uhelp), param.numberOfColumnsIncludingGhost, param.numberOfRowsIncludingGhost, 
													 &gridComm2d, neighborsRanks, &sendWestColumn, &recvWestColumn, &sendEastColumn, &recvEastColumn);
			MPI_Reduce(&localResidual, &globalResidual, 1, MPI_DOUBLE, MPI_SUM, 0, gridComm2d);
		}
		break;

		// omp+mpi hybrid blocking communication
		case 2:
		for (iter = 0; iter < param.maxiter; iter++) {
			localResidual = relax_jacobi_hybrid_blocking(&(param.u), &(param.uhelp), param.numberOfColumnsIncludingGhost, param.numberOfRowsIncludingGhost, 
												  &gridComm2d, neighborsRanks, &sendWestColumn, &recvWestColumn, &sendEastColumn, &recvEastColumn);
			MPI_Reduce(&localResidual, &globalResidual, 1, MPI_DOUBLE, MPI_SUM, 0, gridComm2d);
		}		
		break;

		// omp+mpi hybrid non-blocking communication
		case 3:
		for (iter = 0; iter < param.maxiter; iter++) {
			localResidual = relax_jacobi_hybrid_nonblocking(&(param.u), &(param.uhelp), param.numberOfColumnsIncludingGhost, param.numberOfRowsIncludingGhost, 
													 &gridComm2d, neighborsRanks, &sendWestColumn, &recvWestColumn, &sendEastColumn, &recvEastColumn);
			MPI_Reduce(&localResidual, &globalResidual, 1, MPI_DOUBLE, MPI_SUM, 0, gridComm2d);
		}		
		break;	
		}

		if(rank==0){
			time[exp_number] = wtime() - time[exp_number];
			printf("\n\nResolution: %u\n", param.act_res);
			printf("===================\n");
			printf("Execution time: %f\n", time[exp_number]);
			printf("Residual: %f\n\n", globalResidual);
			printf("megaflops:  %.1lf\n", (double) param.maxiter * param.act_res * param.act_res * 7 / time[exp_number] / 1000000);
			printf("  flop instructions (M):  %.3lf\n", (double) param.maxiter * param.act_res * param.act_res * 7 / 1000000);
		}
		exp_number++;
		free(sendWestColumn);
		free(recvWestColumn);
		free(sendEastColumn);
		free(recvEastColumn);
	}	

	param.act_res = param.act_res - param.res_step_size;
	
	coarsen(param.u, param.act_res + 2, param.act_res + 2, param.startColumnIndex, param.startRowIndex, param.numberOfColumnsIncludingGhost, param.numberOfRowsIncludingGhost, param.uvis, param.visres + 2, param.visres + 2);	

	double * uvisGlobal;
	if(rank==0) uvisGlobal = (double*)calloc( sizeof(double),(param.visres+2)*(param.visres+2) );

	MPI_Reduce(param.uvis, uvisGlobal, (param.visres + 2)*(param.visres + 2), MPI_DOUBLE, MPI_SUM, 0, gridComm2d);
	
	if(rank==0) {

		// for(int i = 0; i < (param.visres + 2); ++i){
		// 	for (int j = 0; j < (param.visres + 2); ++j){
		// 		printf("%f\n", uvisGlobal[i*(param.visres + 2)+j]);
		// 	}
		// }		

		write_image(resfile, uvisGlobal, param.visres + 2, param.visres + 2);
		free(uvisGlobal);
	}
	
	finalize(&param);
	MPI_Finalize();
    return 0;
}