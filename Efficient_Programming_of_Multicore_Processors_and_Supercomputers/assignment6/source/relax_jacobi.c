/*
 * relax_jacobi.c
 *
 * Jacobi Relaxation
 *
 */

#include "heat.h"
#include <omp.h>

double relax_jacobi_blocking( double **u1, double **utmp1,
         unsigned sizex, unsigned sizey,
         MPI_Comm * comm, int * neighborsRanks,
         double** sendWestColumn1, double** recvWestColumn1, 
         double** sendEastColumn1, double** recvEastColumn1 )
{
  int i, j;
  double *help,*u, *utmp,factor=0.5;

  utmp=*utmp1;
  u=*u1;
  double unew, diff, sum=0.0;
  enum directions {north, south, west, east};
  MPI_Status status;
  MPI_Comm gridComm2d = *comm;
  double * sendWestColumn = *sendWestColumn1;
  double * recvWestColumn = *recvWestColumn1;
  double * sendEastColumn = *sendEastColumn1;
  double * recvEastColumn = *recvEastColumn1;

  if (neighborsRanks[south] != MPI_PROC_NULL) {
    MPI_Sendrecv(u + (sizey-2) * sizex + 1, sizex-2, MPI_DOUBLE, neighborsRanks[south], 0,
                 u + (sizey-1) * sizex + 1, sizex-2, MPI_DOUBLE, neighborsRanks[south], 1, gridComm2d, &status);
  }

  if (neighborsRanks[north] != MPI_PROC_NULL) {
    MPI_Sendrecv(u + sizex + 1, sizex-2, MPI_DOUBLE, neighborsRanks[north], 1,
                 u + 1, sizex-2, MPI_DOUBLE, neighborsRanks[north], 0, gridComm2d, &status);
  }

  if (neighborsRanks[east] != MPI_PROC_NULL) {
    for (i = 0; i < sizey-2; i++) {
      sendEastColumn[i] = u[(i + 2) * sizex - 2];
    }
    MPI_Sendrecv(sendEastColumn, sizey-2, MPI_DOUBLE, neighborsRanks[east], 2,
                 recvEastColumn, sizey-2, MPI_DOUBLE, neighborsRanks[east], 3, gridComm2d, &status);
    for (i = 0; i < sizey-2; i++) {
      u[(i + 2) * sizex - 1] = recvEastColumn[i];
    }
  }

  if (neighborsRanks[west] != MPI_PROC_NULL) {
    for (i = 0; i < sizey-2; i++) {
      sendWestColumn[i] = u[(i + 1) * sizex + 1];
    }
    MPI_Sendrecv(sendWestColumn, sizey-2, MPI_DOUBLE, neighborsRanks[west], 3,
                 recvWestColumn, sizey-2, MPI_DOUBLE, neighborsRanks[west], 2, gridComm2d, &status);
    for (i = 0; i < sizey-2; i++) {
      u[(i + 1) * sizex] = recvWestColumn[i];
    }
  }

  for( i=1; i<sizey-1; i++ ) {
  	int ii=i*sizex;
  	int iim1=(i-1)*sizex;
  	int iip1=(i+1)*sizex;
#pragma ivdep
    for( j=1; j<sizex-1; j++ ){
       unew = 0.25 * (u[ ii+(j-1) ]+
        		            u[ ii+(j+1) ]+
        		            u[ iim1+j ]+
        		            u[ iip1+j ]);
		    diff = unew - u[ii + j];
		    utmp[ii+j] = unew;
		    sum += diff * diff;

       }
    }

  *u1=utmp;
  *utmp1=u;
  return(sum);
}

double relax_jacobi_nonblocking( double **u1, double **utmp1,
         unsigned sizex, unsigned sizey,
          MPI_Comm * comm, int * neighborsRanks, 
          double** sendWestColumn1, double** recvWestColumn1, 
          double** sendEastColumn1, double** recvEastColumn1)
{

      int i, j;
      double *help,*u, *utmp,factor=0.5;

      utmp=*utmp1;
      u=*u1;
      double unew, diff, sum=0.0;

      enum directions {north, south, west, east};
      MPI_Request sendRequest[4];
      MPI_Request recvRequest[4];
      MPI_Status status[4];
      MPI_Comm gridComm2d = *comm;
      double * sendWestColumn = *sendWestColumn1;
      double * recvWestColumn = *recvWestColumn1;
      double * sendEastColumn = *sendEastColumn1;
      double * recvEastColumn = *recvEastColumn1;

			if (neighborsRanks[south] != MPI_PROC_NULL){
				MPI_Isend(u + (sizey-2) * sizex + 1, (sizex - 2), MPI_DOUBLE, neighborsRanks[south], 0, gridComm2d, &sendRequest[south]);
				MPI_Irecv(u + (sizey-1) * sizex + 1, (sizex - 2), MPI_DOUBLE, neighborsRanks[south], 1, gridComm2d, &recvRequest[south]);
			}

			if (neighborsRanks[north] != MPI_PROC_NULL){				
				MPI_Isend(u + sizex + 1, (sizex - 2), MPI_DOUBLE, neighborsRanks[north], 1, gridComm2d, &sendRequest[north]);
				MPI_Irecv(u + 1, (sizex - 2), MPI_DOUBLE, neighborsRanks[north], 0, gridComm2d, &recvRequest[north]);
			}

			if (neighborsRanks[east] != MPI_PROC_NULL){
				for (i = 0; i < sizey-2; i++){
					sendEastColumn[i] = u[(i+2)*sizex - 2];
				}
				MPI_Isend(sendEastColumn, (sizey-2), MPI_DOUBLE, neighborsRanks[east], 2, gridComm2d, &sendRequest[east]);
        MPI_Irecv(recvEastColumn, (sizey-2), MPI_DOUBLE, neighborsRanks[east], 3, gridComm2d, &recvRequest[east]);	
      }

			if (neighborsRanks[west] != MPI_PROC_NULL){
				for (i = 0; i < sizey-2; i++){
					sendWestColumn[i] = u[(i+1)*sizex +1];
				}
				MPI_Isend(sendWestColumn, (sizey-2), MPI_DOUBLE, neighborsRanks[west], 3, gridComm2d, &sendRequest[west]);
        MPI_Irecv(recvWestColumn, (sizey-2), MPI_DOUBLE, neighborsRanks[west], 2, gridComm2d, &recvRequest[west]);
      }

  // inner part calculation, no need ghost
  for( i=2; i<sizey-2; i++ ) {
  	int ii=i*sizex;
  	int iim1=(i-1)*sizex;
  	int iip1=(i+1)*sizex;
#pragma ivdep
    for( j=2; j<sizex-2; j++ ){
       unew = 0.25 * (u[ ii+(j-1) ]+
        		            u[ ii+(j+1) ]+
        		            u[ iim1+j ]+
        		            u[ iip1+j ]);
		    diff = unew - u[ii + j];
		    utmp[ii+j] = unew;
		    sum += diff * diff;

       }
    }

  // wait all receive requests
  if (neighborsRanks[east] != MPI_PROC_NULL){
    MPI_Wait(&recvRequest[east], &status[east]); 
  }
  if (neighborsRanks[west] != MPI_PROC_NULL){
    MPI_Wait(&recvRequest[west], &status[west]); 
  }
  if (neighborsRanks[south] != MPI_PROC_NULL){
    MPI_Wait(&recvRequest[south], &status[south]); 
  }   
  if (neighborsRanks[north] != MPI_PROC_NULL){
    MPI_Wait(&recvRequest[north], &status[north]); 
  }

  // east memory copy to ghost
  if (neighborsRanks[east] != MPI_PROC_NULL){
    for (i = 0; i < sizey-2; i++){
      u[(i+2)*sizex - 1] = recvEastColumn[i];
    }
  }

  // west memory copy to ghost
  if (neighborsRanks[west] != MPI_PROC_NULL){
    for (i = 0; i < sizey - 2; i++){
      u[(i+1)*sizex] = recvWestColumn[i];
    }
  }

  // sub-east column (not ghost)
  #pragma ivdep
  for (i = 2; i < sizey - 2; i++){
      int ii = i * sizex;
      int iim1 = (i - 1) * sizex;
      int iip1 = (i + 1) * sizex;
      j = sizex-2;
        unew = 0.25 * (u[ii + (j - 1)] +
                        u[ii + (j + 1)] +
                        u[iim1 + j] +
                        u[iip1 + j]);
        diff = unew - u[ii + j];
        utmp[ii + j] = unew;
        sum += diff * diff;
    }

  // sub-west column (not ghost)
  #pragma ivdep
  for (i = 2; i < sizey - 2; i++){
    int ii = i * sizex;
    int iim1 = (i - 1) * sizex;
    int iip1 = (i + 1) * sizex;

      unew = 0.25 * (u[ii] +
                      u[ii + 2] +
                      u[iim1+1] +
                      u[iip1+1]);
      diff = unew - u[ii+1];
      utmp[ii + 1] = unew;
      sum += diff * diff;
  }

  // sub-north row (not ghost)
  i = 1;
  int ii=i*sizex;
  int iim1=(i-1)*sizex;
  int iip1=(i+1)*sizex;   
  #pragma ivdep
  for (j = 1; j < sizex - 1; j++){
    unew = 0.25 * (u[ii + (j - 1)] +
                    u[ii + (j + 1)] +
                    u[iim1 + j] +
                    u[iip1 + j]);
    diff = unew - u[ii + j];
    utmp[ii + j] = unew;
    sum += diff * diff;
  }

  // sub-south row (not ghost)
  i = sizey - 2;
  ii=i*sizex;
  iim1=(i-1)*sizex;
  iip1=(i+1)*sizex; 
  #pragma ivdep
  for (j = 1; j < sizex - 1; j++){
    unew = 0.25 * (u[ii + (j - 1)] +
                    u[ii + (j + 1)] +
                    u[iim1 + j] +
                    u[iip1 + j]);
    diff = unew - u[ii + j];
    utmp[ii + j] = unew;
    sum += diff * diff;
  }

  *u1=utmp;
  *utmp1=u;
  return(sum);
}

double relax_jacobi_hybrid_blocking( double **u1, double **utmp1,
         unsigned sizex, unsigned sizey,
         MPI_Comm * comm, int * neighborsRanks,
         double** sendWestColumn1, double** recvWestColumn1, 
         double** sendEastColumn1, double** recvEastColumn1 )
{
  int i, j;
  double *help,*u, *utmp,factor=0.5;

  utmp=*utmp1;
  u=*u1;
  double unew, diff, sum=0.0;
  enum directions {north, south, west, east};
  MPI_Status status;
  MPI_Comm gridComm2d = *comm;
  double * sendWestColumn = *sendWestColumn1;
  double * recvWestColumn = *recvWestColumn1;
  double * sendEastColumn = *sendEastColumn1;
  double * recvEastColumn = *recvEastColumn1;

  if (neighborsRanks[south] != MPI_PROC_NULL) {
    MPI_Sendrecv(u + (sizey-2) * sizex + 1, sizex-2, MPI_DOUBLE, neighborsRanks[south], 0,
                 u + (sizey-1) * sizex + 1, sizex-2, MPI_DOUBLE, neighborsRanks[south], 1, gridComm2d, &status);
  }

  if (neighborsRanks[north] != MPI_PROC_NULL) {
    MPI_Sendrecv(u + sizex + 1, sizex-2, MPI_DOUBLE, neighborsRanks[north], 1,
                 u + 1, sizex-2, MPI_DOUBLE, neighborsRanks[north], 0, gridComm2d, &status);
  }

  if (neighborsRanks[east] != MPI_PROC_NULL) {
    for (i = 0; i < sizey-2; i++) {
      sendEastColumn[i] = u[(i + 2) * sizex - 2];
    }
    MPI_Sendrecv(sendEastColumn, sizey-2, MPI_DOUBLE, neighborsRanks[east], 2,
                 recvEastColumn, sizey-2, MPI_DOUBLE, neighborsRanks[east], 3, gridComm2d, &status);
    for (i = 0; i < sizey-2; i++) {
      u[(i + 2) * sizex - 1] = recvEastColumn[i];
    }
  }

  if (neighborsRanks[west] != MPI_PROC_NULL) {
    for (i = 0; i < sizey-2; i++) {
      sendWestColumn[i] = u[(i + 1) * sizex + 1];
    }
    MPI_Sendrecv(sendWestColumn, sizey-2, MPI_DOUBLE, neighborsRanks[west], 3,
                 recvWestColumn, sizey-2, MPI_DOUBLE, neighborsRanks[west], 2, gridComm2d, &status);
    for (i = 0; i < sizey-2; i++) {
      u[(i + 1) * sizex] = recvWestColumn[i];
    }
  }
#pragma omp parallel for private(j, unew, diff) reduction(+:sum) schedule(static)
  for( i=1; i<sizey-1; i++ ) {
  	int ii=i*sizex;
  	int iim1=(i-1)*sizex;
  	int iip1=(i+1)*sizex;
#pragma ivdep
    for( j=1; j<sizex-1; j++ ){
       unew = 0.25 * (u[ ii+(j-1) ]+
        		            u[ ii+(j+1) ]+
        		            u[ iim1+j ]+
        		            u[ iip1+j ]);
		    diff = unew - u[ii + j];
		    utmp[ii+j] = unew;
		    sum += diff * diff;

       }
    }

  *u1=utmp;
  *utmp1=u;
  return(sum);
}

double relax_jacobi_hybrid_nonblocking( double **u1, double **utmp1,
         unsigned sizex, unsigned sizey, 
         MPI_Comm * comm, int * neighborsRanks, 
          double** sendWestColumn1, double** recvWestColumn1, 
          double** sendEastColumn1, double** recvEastColumn1)
{

      int i, j;
      double *help,*u, *utmp,factor=0.5;

      utmp=*utmp1;
      u=*u1;
      double unew, diff, sum=0.0;

      enum directions {north, south, west, east};
      MPI_Request sendRequest[4];
      MPI_Request recvRequest[4];
      MPI_Status status[4];
      MPI_Comm gridComm2d = *comm;
      double * sendWestColumn = *sendWestColumn1;
      double * recvWestColumn = *recvWestColumn1;
      double * sendEastColumn = *sendEastColumn1;
      double * recvEastColumn = *recvEastColumn1;

			if (neighborsRanks[south] != MPI_PROC_NULL){
				MPI_Isend(u + (sizey-2) * sizex + 1, (sizex - 2), MPI_DOUBLE, neighborsRanks[south], 0, gridComm2d, &sendRequest[south]);
				MPI_Irecv(u + (sizey-1) * sizex + 1, (sizex - 2), MPI_DOUBLE, neighborsRanks[south], 1, gridComm2d, &recvRequest[south]);
			}

			if (neighborsRanks[north] != MPI_PROC_NULL){				
				MPI_Isend(u + sizex + 1, (sizex - 2), MPI_DOUBLE, neighborsRanks[north], 1, gridComm2d, &sendRequest[north]);
				MPI_Irecv(u + 1, (sizex - 2), MPI_DOUBLE, neighborsRanks[north], 0, gridComm2d, &recvRequest[north]);
			}

			if (neighborsRanks[east] != MPI_PROC_NULL){
				for (i = 0; i < sizey-2; i++){
					sendEastColumn[i] = u[(i+2)*sizex - 2];
				}
				MPI_Isend(sendEastColumn, (sizey-2), MPI_DOUBLE, neighborsRanks[east], 2, gridComm2d, &sendRequest[east]);
        MPI_Irecv(recvEastColumn, (sizey-2), MPI_DOUBLE, neighborsRanks[east], 3, gridComm2d, &recvRequest[east]);	
      }

			if (neighborsRanks[west] != MPI_PROC_NULL){
				for (i = 0; i < sizey-2; i++){
					sendWestColumn[i] = u[(i+1)*sizex +1];
				}
				MPI_Isend(sendWestColumn, (sizey-2), MPI_DOUBLE, neighborsRanks[west], 3, gridComm2d, &sendRequest[west]);
        MPI_Irecv(recvWestColumn, (sizey-2), MPI_DOUBLE, neighborsRanks[west], 2, gridComm2d, &recvRequest[west]);
      }

  // inner part calculation, no need ghost
#pragma omp parallel for private(j, unew, diff) reduction(+:sum) schedule(static)
  for( i=2; i<sizey-2; i++ ) {
  	int ii=i*sizex;
  	int iim1=(i-1)*sizex;
  	int iip1=(i+1)*sizex;
#pragma ivdep
    for( j=2; j<sizex-2; j++ ){
       unew = 0.25 * (u[ ii+(j-1) ]+
        		            u[ ii+(j+1) ]+
        		            u[ iim1+j ]+
        		            u[ iip1+j ]);
		    diff = unew - u[ii + j];
		    utmp[ii+j] = unew;
		    sum += diff * diff;

       }
    }

  // wait all receive requests
  if (neighborsRanks[east] != MPI_PROC_NULL){
    MPI_Wait(&recvRequest[east], &status[east]); 
  }
  if (neighborsRanks[west] != MPI_PROC_NULL){
    MPI_Wait(&recvRequest[west], &status[west]); 
  }
  if (neighborsRanks[south] != MPI_PROC_NULL){
    MPI_Wait(&recvRequest[south], &status[south]); 
  }   
  if (neighborsRanks[north] != MPI_PROC_NULL){
    MPI_Wait(&recvRequest[north], &status[north]); 
  }

  // east memory copy to ghost
  if (neighborsRanks[east] != MPI_PROC_NULL){
    for (i = 0; i < sizey-2; i++){
      u[(i+2)*sizex - 1] = recvEastColumn[i];
    }
  }

  // west memory copy to ghost
  if (neighborsRanks[west] != MPI_PROC_NULL){
    for (i = 0; i < sizey - 2; i++){
      u[(i+1)*sizex] = recvWestColumn[i];
    }
  }

  // sub-east column (not ghost)
#pragma omp parallel for private(j, unew, diff) reduction(+:sum) schedule(static)
  for (i = 2; i < sizey - 2; i++){
      int ii = i * sizex;
      int iim1 = (i - 1) * sizex;
      int iip1 = (i + 1) * sizex;
      j = sizex-2;
        unew = 0.25 * (u[ii + (j - 1)] +
                        u[ii + (j + 1)] +
                        u[iim1 + j] +
                        u[iip1 + j]);
        diff = unew - u[ii + j];
        utmp[ii + j] = unew;
        sum += diff * diff;
    }

  // sub-west column (not ghost)
#pragma omp parallel for private(unew, diff) reduction(+:sum) schedule(static)
  for (i = 2; i < sizey - 2; i++){
    int ii = i * sizex;
    int iim1 = (i - 1) * sizex;
    int iip1 = (i + 1) * sizex;

      unew = 0.25 * (u[ii] +
                      u[ii + 2] +
                      u[iim1+1] +
                      u[iip1+1]);
      diff = unew - u[ii+1];
      utmp[ii + 1] = unew;
      sum += diff * diff;
  }

  // sub-north row (not ghost)
  i = 1;
  int ii=i*sizex;
  int iim1=(i-1)*sizex;
  int iip1=(i+1)*sizex;   
#pragma omp parallel for private(unew, diff) reduction(+:sum) schedule(static)
  for (j = 1; j < sizex - 1; j++){
    unew = 0.25 * (u[ii + (j - 1)] +
                    u[ii + (j + 1)] +
                    u[iim1 + j] +
                    u[iip1 + j]);
    diff = unew - u[ii + j];
    utmp[ii + j] = unew;
    sum += diff * diff;
  }

  // sub-south row (not ghost)
  i = sizey - 2;
  ii=i*sizex;
  iim1=(i-1)*sizex;
  iip1=(i+1)*sizex; 
#pragma omp parallel for private(unew, diff) reduction(+:sum) schedule(static)
  for (j = 1; j < sizex - 1; j++){
    unew = 0.25 * (u[ii + (j - 1)] +
                    u[ii + (j + 1)] +
                    u[iim1 + j] +
                    u[iip1 + j]);
    diff = unew - u[ii + j];
    utmp[ii + j] = unew;
    sum += diff * diff;
  }

  *u1=utmp;
  *utmp1=u;
  return(sum);
}