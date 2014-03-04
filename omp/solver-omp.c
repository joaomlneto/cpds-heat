#include "heat.h"

#define NB 8

#define min(a,b) ( ((a) < (b)) ? (a) : (b) )

/*
 * Blocked Jacobi solver: one iteration step
 */
double relax_jacobi (double *u, double *utmp, unsigned sizex, unsigned sizey)
{
	double diff, sum=0.0;
	int nbx, bx, nby, by;

	nbx = omp_get_max_threads();//NB;
	bx = sizex/nbx + ((sizex%nbx) ? 1 : 0);//sizex/nbx;
	nby = 1;//NB;
	by = sizey/nby;

	#pragma omp parallel for reduction(+:sum) private(diff)
	for (int ii=0; ii<nbx; ii++) {
		for (int jj=0; jj<nby; jj++)  {
			for (int i=1+ii*bx; i<=min((ii+1)*bx, sizex-2); i++) {
				for (int j=1+jj*by; j<=min((jj+1)*by, sizey-2); j++) {
					utmp[i*sizey+j]= 0.25 * (u[ i*sizey	    + (j-1) ]+  // left
					                         u[ i*sizey     + (j+1) ]+  // right
					                         u[ (i-1)*sizey + j     ]+  // top
					                         u[ (i+1)*sizey + j     ]); // bottom
					diff = utmp[i*sizey+j] - u[i*sizey + j];
					sum += diff * diff; 
				}
			}
		}
	}

	return sum;
}

/*
 * Blocked Red-Black solver: one iteration step
 */
double relax_redblack (double *u, unsigned sizex, unsigned sizey)
{
	double unew, diff, sum=0.0;
	int nbx, bx, nby, by;
	int lsw;

//*
	nbx = NB;
	bx = sizex/nbx;
	nby = NB;
	by = sizey/nby;
// */

/*
	nbx = omp_get_max_threads();
	bx = sizex/nbx + ((sizex%nbx) ? 1 : 0);
	nby = 1;
	by = sizey/nby;
// */

	// Computing "Red" blocks
	#pragma omp parallel for reduction(+:sum) private(diff, lsw)
	for (int ii=0; ii<nbx; ii++) {
		lsw = ii%2;
		for (int jj=lsw; jj<nby; jj=jj+2) {
			for (int i=1+ii*bx; i<=min((ii+1)*bx, sizex-2); i++) {
				for (int j=1+jj*by; j<=min((jj+1)*by, sizey-2); j++) {
					unew= 0.25 * (	u[ i*sizey     + (j-1) ]+  // left
					               u[ i*sizey     + (j+1) ]+  // right
					               u[ (i-1)*sizey	+ j     ]+  // top
					               u[ (i+1)*sizey	+ j     ]); // bottom
					diff = unew - u[i*sizey+ j];
					sum += diff * diff;
					u[i*sizey+j]=unew;
				}
			}
		}
	}

	#pragma omp parallel for reduction(+:sum) private(diff, lsw)
	// Computing "Black" blocks
	for (int ii=0; ii<nbx; ii++) {
		lsw = (ii+1)%2;
		for (int jj=lsw; jj<nby; jj=jj+2) {
			for (int i=1+ii*bx; i<=min((ii+1)*bx, sizex-2); i++) {
				for (int j=1+jj*by; j<=min((jj+1)*by, sizey-2); j++) {
					unew= 0.25 * (	u[ i*sizey     + (j-1) ]+  // left
					               u[ i*sizey     + (j+1) ]+  // right
					               u[ (i-1)*sizey	+ j     ]+  // top
					               u[ (i+1)*sizey	+ j     ]); // bottom
					diff = unew - u[i*sizey+ j];
					sum += diff * diff; 
					u[i*sizey+j]=unew;
				}
			}
		}
	}

	return sum;
}

/*
 * Blocked Gauss-Seidel solver: one iteration step
 */
double relax_gauss (double *u, unsigned sizex, unsigned sizey)
{
	double unew, diff, sum=0.0;
	int nbx, bx, nby, by;

/*
	nbx = NB;
	bx = sizex/nbx;
	nby = NB;
	by = sizey/nby;
// */

//*
	nbx = omp_get_max_threads();
	bx = sizex/nbx + ((sizex%nbx) ? 1 : 0);
	nby = 1;
	by = sizey/nby;
// */

	#pragma omp parallel for reduction(+:sum) private(diff)
	for (int ii=0; ii<nbx; ii++)
		for (int jj=0; jj<nby; jj++) 
			for (int i=1+ii*bx; i<=min((ii+1)*bx, sizex-2); i++) 
				for (int j=1+jj*by; j<=min((jj+1)*by, sizey-2); j++) {
					unew= 0.25 * (	u[ i*sizey     + (j-1) ]+  // left
					               u[ i*sizey     + (j+1) ]+  // right
					               u[ (i-1)*sizey	+ j     ]+  // top
					               u[ (i+1)*sizey	+ j     ]); // bottom
					diff = unew - u[i*sizey+ j];
					sum += diff * diff; 
					u[i*sizey+j]=unew;
				}

	return sum;
}

