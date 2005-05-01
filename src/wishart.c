#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include "latentUtil.h"

void WishartRandom (long n, long p, double *t, double *w, double *work);
void rwish_wrapper (long *n, long *p, double *t, double *w);
void riwish_wrapper (long *n, long *p, double *t, double *w, double **dd1,
		     double **dd2, double *vdd2);

void call_riwish (long n, long p, double *t, double **w, double **dd1,
		  double **dd2, double *vdd, double *vdd2)
{
  double *mat1;
  int i,j;

  /*   mat1 = dvector(p * p); */
  mat1 = vdd;

  /*   Rprintf("In riwish\n"); */

  for(i=0;i<p;i++)
    for(j=0;j<p;j++)
      mat1[i +j *p] = w[i][j];

  /*   Rprintf("riwish_wrapper\n"); */
  riwish_wrapper (&n, &p, t, mat1, dd1, dd2,vdd2);
  /*   Rprintf("OUT riwish_wrapper\n"); */

  for(i=0;i<p;i++)
    for(j=0;j<p;j++)
      w[i][j] = mat1[i +j *p];

  /*   free_dvector(mat1); */
}


void rwish_wrapper (long *n, long *p, double *t, double *w)
{
  double *work;

  work = (double *) malloc(*p * *p * sizeof(double));
  GetRNGstate(); 
  WishartRandom (*n, *p, t, w, work);
  PutRNGstate();
  free(work);
}

void riwish_wrapper (long *n, long *p, double *t, double *w, double **dd1,
		     double **dd2, double *vdd2)
{
  double *work, **mat1, **mat2, sum;
  int i,j,k,pint;
  /*   GetRNGstate();  */

  /*   work = (double *) malloc(*p * *p * sizeof(double)); */
  work = vdd2;

  /*   Rprintf("In riwish_wrapper\n"); */
  pint = (int)*p;

  /*   mat1 = dmatrix(*p,*p); */
  /*   mat2 = dmatrix(*p,*p); */
  mat1 = dd1;
  mat2 = dd2;

  /*   Rprintf("WishartRandom\n"); */
  WishartRandom (*n, *p, t, w, work);
  /*   Rprintf("OUT WishartRandom\n"); */
  k = 0;
  for(j=0; j< *p; j++)
    for(i=0; i<= j; i++)
    {
      mat1[i][j]= mat1[j][i] = w[k];
      k++;
    }
  Rprintf("blah = [%1.4g, %1.4g; %1.4g, %1.4g]\n",mat1[0][0],mat1[0][1],mat1[1][0],mat1[1][1]);
  /*   Rprintf("Inverse!!!!!!!!!!!!!!\n"); */
  if(*p==2)
  {
    sum = mat1[0][0] * mat1[1][1] - mat1[0][1] * mat1[1][0];
    mat2[0][0] = mat1[1][1]/sum;
    mat2[1][1] = mat1[0][0]/sum;
    mat2[1][0] = mat2[0][1] = -mat1[1][0]/sum;
  }
  else
    inverse(mat1,&pint,mat2);
  /*   Rprintf("OUT Inverse\n"); */
  Rprintf("blah = [%1.4g, %1.4g; %1.4g, %1.4g]\n",mat2[0][0],mat2[0][1],mat2[1][0],mat2[1][1]);

  k = 0;
  for(j=0; j< *p; j++)
    for(i=0; i< *p; i++)
    {
      w[k] = mat2[i][j];
      k++;
    }
  /*   PutRNGstate(); */
  /*   free(work); free_dmatrix(mat1,*p); free_dmatrix(mat2,*p); */
}


/*
-----------------------------------------------------------------------------
   Wishart

   samples a matrix according to the Wishart distribution by the method
   of Odell and Feiveson (1966).

   Parameters are:
   n (degrees of freedom); p (dimension of Wishart matrix);
   t (pointer to a Cholesky decomposition of a covariance matrix);
   w (pointer to the sampled Wishart matrix, in
   triangular form; work (pointer to a work space, of length p*p).

   Triangular matrices are stored in order
   0 1 3
     2 4
       5 etc.
*/

void WishartRandom (long n, long p, double *t, double *w, double *work)
{
  double eta, sum;
  long i, j, k, m, k1, k2, k3;

  /*   Rprintf ("::WishartRandom::\nn = %d :: p = %d\n",n,p); */
  /*   exit(0); */

  /* generate random quantities for Bartlett's decomposition */
  for (j = 0, k = 0; j < p; j++) {
    for (i = 0; i < j; i++)
      w[k++] = rnorm(0, 1);
    /* Chi-square with n-i degrees of freedom */
    eta = ( ((double)n) - ((double)i) ) / 2.0;
    /*     Rprintf("j = %d || i = %d || k = %d || rnchisq(%1.4f,0.5)\n",j,i,k,eta); */
    sum = rnchisq(eta, 0.5);
    /*     Rprintf("done the chisq\n"); */
    w[k++] = sum;
    /*     Rprintf("done the chisqthingee\n"); */
  }
  /*   Rprintf("1--\n"); */

  /* generate a standard Wishart */
  for (j = p - 1, m = k - 1, k2 = (p * (p - 1)) / 2; j >= 0; k2 = k2 - (j--)) {
    eta = w[m];
    for (i = j, k1 = (i * (i + 1)) / 2; i >= 0; k1 = k1 - (i--), m--) {
      for (k = 0, sum = 0.0; k < i; k++)
        sum = sum + w[k1+k] * w[k2+k];

      if (i == j)
        w[m] = sum + eta;
      else
        w[m] = sum + sqrt(eta) * w[m];
    }
  }

  /*   Rprintf("2--\n"); */

  /* form product L * W * L' */
  for (i = 0, k1 = 0, m = 0; i < p; k1 = k1 + (++i)) {
    for (j = 0, k2 = 0; j < p; k2 = k2 + (++j), m++) {
      for (k = 0, sum = 0.0; k < j; k++)
        sum = sum + t[k1+k] * w[k2+k];

      for (k = j, k3 = j; k <= i; k3 = k3 + (++k))
        sum = sum + t[k1+k] * w[k2+k3];
     for (k = j, k3 = j; k <= i; k3 = k3 + (++k))
        sum = sum + t[k1+k] * w[k2+k3];

      work[m] = sum;
    }
  }

  /*   Rprintf("3--\n"); */

  for (i = 0, m = 0, k1 = 0; i < p; i++, k1 = k1 + p) {
    for (j = 0, k2 = 0; j <= i; k2 = k2 + (++j), m++) {
      for (k = 0, sum = 0.0; k <= j; k++)
        sum = sum + work[k1+k] * t[k2+k];

      w[m] = sum;
    }
  }
  /*   Rprintf("OUT-->\n"); */

} /* WishartRandom */



/* Dirichlet generator */
void rdirich(int n, double *epsilon) 
{
  int i;
  double *y, ysum = 0;
  /*   Rprintf("n = %d\n",n); */
  y = dvector(n);
  for (i=0; i<n; ++i){
    /*     if(epsilon[i]==1)epsilon[i] +=0.5; */
    /*     Rprintf("epsilon[%d] = %1.7f\n",i,epsilon[i]); */
    y[i] = rgamma(epsilon[i], 1.0);
    ysum += y[i];
  }
  
  for (i=0; i<n; ++i)
    epsilon[i] = y[i]/ysum;
  free(y);
}

void dirichlet_wrapper(int *n, double *epsilon)
{
  GetRNGstate(); 
  rdirich(*n, epsilon);
  PutRNGstate();
}
