/*
 * relax_jacobi.c
 *
 * Jacobi Relaxation
 *
 */

#include "heat.h"

/*
 * One Jacobi iteration step including residual caculation
 */
double relax_jacobi(double *u, double *utmp, unsigned sizex, unsigned sizey) {
	unsigned i, j;
	double diff, sum = 0.0;
	for (i = 1; i < sizey - 1; i++) {
		for (j = 1; j < sizex - 1; j++) {
			utmp[i * sizex + j] = 0.25 * (u[i * sizex + (j - 1)] +  // left
						u[i * sizex + (j + 1)] +  // right
						u[(i - 1) * sizex + j] +  // top
						u[(i + 1) * sizex + j]); // bottom
						diff = utmp[i * sizex + j] - u[i * sizex + j];
						sum += diff * diff;
		}
	}
	return sum;
}