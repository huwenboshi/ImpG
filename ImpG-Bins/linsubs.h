#include  <stdio.h>
#include <string.h>
#include <unistd.h>
#include <math.h>

#define MAX(x,y) ((x)>(y)?(x):(y))
#define SIGN(a, b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

#define ROTATE(a,i,j,k,l) g=a[i*n+j];h=a[k*(n)+l];a[i*(n)+j]=g-s*(h+g*tau);\
	a[k*(n)+l]=h+s*(g-h*tau);

/* linear algebra */
void pdinv(double *cinv, double *coeff, int n) ;

// compute pseudo inverse
void pinv(double *A_inv, double *A, size_t n);

// compute pseudo inverse using jacobi
void pinv_jacobi(double *A_inv, double *A, size_t n);

/* numer recipes p 97 */
int choldc (double *a, int n, double p[]);
void cholsl (double *a, int n, double p[], double b[], double x[]);
void cholesky(double *cf, double *a, int n) ;
void pmat(double *mat, int n)   ;
