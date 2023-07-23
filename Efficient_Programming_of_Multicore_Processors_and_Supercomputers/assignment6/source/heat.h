/*
 * heat.h
 *
 * Global definitions for the iterative solver
 */

#ifndef JACOBI_H_INCLUDED
#define JACOBI_H_INCLUDED

#include <stdio.h>
#include <mpi.h>
// configuration

typedef struct
{
    float posx;
    float posy;
    float range;
    float temp;
}
heatsrc_t;

typedef struct
{
    unsigned maxiter;       // maximum number of iterations
    unsigned act_res;
    unsigned max_res;       // spatial resolution
    unsigned initial_res;
    unsigned res_step_size;
    unsigned visres;        // visualization resolution
  
    double *u, *uhelp;
    double *uvis;
    
    int algorithm;

    int numberOfRows;
    int numberOfColumns;
    int startRowIndex;
    int startColumnIndex;
    int endRowIndex;
    int endColumnIndex;
    int numberOfRowsIncludingGhost;
    int numberOfColumnsIncludingGhost;
    unsigned   numsrcs;     // number of heat sources
    heatsrc_t *heatsrcs;
}
algoparam_t;


// function declarations

// misc.c
int initialize( algoparam_t *param );
int finalize( algoparam_t *param );
void write_image( FILE * f, double *u,
		  unsigned sizex, unsigned sizey );
int coarsen(double *uold, unsigned oldx, unsigned oldy , 
        int startColumnIndex, int startRowIndex, int columns, int rows,
	    double *unew, unsigned newx, unsigned newy );

// Gauss-Seidel: relax_gauss.c
double residual_gauss( double *u, double *utmp,
		       unsigned sizex, unsigned sizey );
void relax_gauss( double *u, 
		  unsigned sizex, unsigned sizey  );

// Jacobi: relax_jacobi.c
double residual_jacobi( double *u,
			unsigned sizex, unsigned sizey );
double relax_jacobi_blocking( double **u1, double **utmp1,
         unsigned sizex, unsigned sizey,
         MPI_Comm * comm, int * neighborsRanks,
         double** sendWestColumn1, double** recvWestColumn1, 
         double** sendEastColumn1, double** recvEastColumn1 );
double relax_jacobi_nonblocking( double **u1, double **utmp1,
         unsigned sizex, unsigned sizey, 
         MPI_Comm * comm, int * neighborsRanks, 
          double** sendWestColumn1, double** recvWestColumn1, 
          double** sendEastColumn1, double** recvEastColumn1); 
double relax_jacobi_hybrid_blocking( double **u1, double **utmp1,
         unsigned sizex, unsigned sizey,
         MPI_Comm * comm, int * neighborsRanks,
         double** sendWestColumn1, double** recvWestColumn1, 
         double** sendEastColumn1, double** recvEastColumn1 );
double relax_jacobi_hybrid_nonblocking( double **u1, double **utmp1,
         unsigned sizex, unsigned sizey, 
         MPI_Comm * comm, int * neighborsRanks, 
          double** sendWestColumn1, double** recvWestColumn1, 
          double** sendEastColumn1, double** recvEastColumn1);          
#endif // JACOBI_H_INCLUDED
