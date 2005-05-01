/****************************************************************************/
/*  Author: Susan Shortreed, susanms@stat.washington.edu                    */
/*  Purpose: support functions for parameter estimation for the latent      */
/*           space model proposed by  Peter D. Hoff, Adrian E. Raftery      */
/*           and Mark S. Handcock in                                        */
/*           "Latent Space Approaches to Social Network Analysis"           */
/*           All of this code is for an R function to be incorporated       */
/*           into the R ERGMM package.                                       */
/****************************************************************************/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <R.h>
#include <Rmath.h>
#include "latentUtil.h"

#define FOUND 0
#define NOTFOUND 1



/*  General Notes
    n: number of rows in a matrix or length of a vector
    m: number of colums in a matrix
*/

/*  Allocates memory for a vector of doubles of length n */
double *dvector(int n){
  double *a;
  int i;
  a = (double*) malloc(n*sizeof(double));
  if(a == NULL) Rprintf("Not enough memory to make double vector");
  for(i=0;i<n;i++){
    a[i]=0.0;
  }
  return a;
}

/*  Allocates memory for a vector of doubles of length n */
int *ivector(int n){
  int *a;
  int i;
  a = (int*) malloc(n*sizeof(int));
  if(a == NULL) Rprintf("Not enough memory to make integer vector");
  for(i=0;i<n;i++){
    a[i]=0;
  }
  return a;
}

/*  Allocates memory for an n by m matrix of doubles */
double **dmatrix(int n,int m)
{
  double **A;
  int i, j;

  /* assigning memory and initialize */
  A = (double**) malloc(n*sizeof(double*));
  if(A == NULL) Rprintf("Not enough memory to make double matrix");
  for(i=0;i<n;i++){
     A[i] = (double*) malloc(m*sizeof(double));
     if(A[i] == NULL) Rprintf("Not enough memory to make double matrix");
     for(j=0;j<m;j++){
       A[i][j]=0.0;
     }
  }

  return A;
}


/*  Allocates memory for an n by m matrix of doubles */
int **imatrix(int n,int m)
{
  int **A;
  int i, j;
  /* assigning memory and initialize */
  A = (int**) malloc(n*sizeof(int*));
  if(A == NULL) Rprintf("Not enough memory to make integer matrix");
  for(i=0;i<n;i++){
     A[i] = (int*) malloc(m*sizeof(int));
     if(A[i] == NULL) Rprintf("Not enough memory to make integer matrix");
     for(j=0;j<m;j++){
       A[i][j]=0;
     }
  }
  return A;
}

/*  Frees memory for a vector of doubles of length n */
void free_dvector(double *x){
    free(x);
}

/*  Frees memory for a vector of doubles of length n */
void free_ivector(int *x){
    free(x);
}

/* Frees memory for an n by m matrix of doubles*/
void free_dmatrix(double **A, int n){
  int i=0;
  for(i=0;i<n;i++){
    free(A[i]);
  }
  free(A);
}

/* Frees memory for an n by m matrix of doubles*/
void free_imatrix(int **A, int n){
  int i=0;
  for(i=0;i<n;i++){
    free(A[i]);
  }
  free(A);
}

void init_dvector(double *x, int n, double value)
{
  int i;
  for(i=0;i<n;i++){
    x[i] = value;
  }
}

void init_dmatrix(double **A, int n, int m, double value)
{
  int i,j;
  for(i=0;i<n;i++){
    for(j=0;j<m;j++){
      A[i][j]=value;
    }
  }
}


/* Concatenates two vectors, the second on the bottom */
double *cat_dvectors(double *x, int nx, double *y, int ny)
{
  int i;
  double *xy = dvector(nx+ny);

  for(i=0;i<nx;i++){
    xy[i] = x[i];
  }
  for(i=0;i<ny;i++){
    xy[i+nx] = y[i];
  }
  return(xy);
}

/* Concatenates a vector and a scaler */
/* if end equals true then place at end */
double *cat_dvector_scaler(double *x, int nx, double y, int end)
{
  int i;
  double *xy; 
  xy = dvector(nx+1);

  if(end==0){
    for(i=0;i<nx;i++){
      xy[i] = x[i];
    }
    xy[nx] = y;
  }
  else{
    xy[0] = y;
    for(i=1;i<=nx;i++){
      xy[i] = x[i-1];
    }
  }
  return(xy);
}


/* Computes AX, where A is na by ma and B is nb by mb */
/*  b needs ot be a vector of length m*/
double *dvector_times_matrix(double *x, int n, double **A, int m, double *b)
{
  
  int i=0,j=0;
  for(i=0;i<m;i++){
    for(j=0;j<n;j++){
      b[i] += x[j]*A[j][i];
    }
  }
  return(b);
}

/* mutliply the scaler times the matrix and return the result in B*/
void dscaler_times_matrix(double x, double **A, int n, int m, double **B)
{
  int i,j;
  for(i=0;i<n;i++){
    for(j=0;j<m;j++){
      B[i][j] = A[i][j] * x;
    } 
  }
}

/* Computes AB, where A is na by ma and B is ma by mb */
void dmatrix_multiply(double **A,int na,int ma, double **B, int mb, 
			  double **C)
{
  int i=0,j=0,k=0;
 
  for(i=0;i<na;i++){
    for(j=0;j<mb;j++){
      for(k=0;k<ma;k++){
	C[i][j] += A[i][k]*B[k][j];
      }
    }
  }
  return;
}

/* Computes AB, where A is na by ma and B is ma by mb */
void imatrix_multiply(int **A,int na,int ma, int **B, int mb, int **C)
{
  int i=0,j=0,k=0;
 
  /* if(ma != nb) return (double**)(-1); */
  for(i=0;i<na;i++){
    for(j=0;j<mb;j++){
      for(k=0;k<ma;k++){
	C[i][j] += A[i][k]*B[k][j];
      }
    }
  }
  return;
}

/* add the second matrix to the first */
void dmatrix_addition(double **A, int n, int m, double **B)
{
  int i,j;
  for(i=0;i<n;i++){
    for(j=0;j<m;j++){
      A[i][j] += B[i][j];
    }
  }
}


/* Makes sure that all values of the matrix are the same*/
void init_imatrix(int **A, int n, int m, int value)
{
  int i,j;

  for(i=0;i<n;i++)
    for(j=0;j<m;j++)
      A[i][j] = value;
  return;
}


/* Returns the transpose of the matrix */
void t(double **A, int n, int m, double **tA)
{
   int i,j;
   for(i=0;i<n;i++){
     for(j=0;j<m;j++){
       tA[j][i] = A[i][j];
     }  
   }
   return;
}


void copy_dmatrix(double **source,double **dest,int n,int m)
{
  int i,j;
  for(i=0;i<n;i++)
    for(j=0;j<m;j++)
      dest[i][j]=source[i][j];
  return;
}

double mean(double *x, int n){
  int i=0;
  double mu=0.0;
  for(i=0;i<n;i++){
    mu += x[i];
  }
  return(mu/n);
}

/* Thanks to Rahpael Gottardo*/
void qr_solve(double **x, int *n1, double ** y, double **coef, int singular)
/* Translation of the R function qr.solve into pure C
   NB We have to transpose the matrices since the ordering of an array is different in Fortran
   NB2 We have to copy x to avoid it being overwritten.
*/

{
    int i,j, info = 0, rank, *pivot, n, p;
    double tol = 1.0E-7, *qraux, *work;
    double * xt, *yt, *coeft;

    qraux = dvector(*n1);
    pivot = ivector(*n1);
    work  = dvector(2*(*n1));
    
    for(i = 0; i < *n1; i++)
        pivot[i] = i+1;

    /** Copy the matrix by column **/
    xt = dvector((*n1)*(*n1));
    for(i=0;i<*n1;i++)
      for(j=0;j<*n1;j++)
	xt[i*(*n1)+j]=x[j][i];

    n = *n1;
    p = *n1;
    
    F77_CALL(dqrdc2)(xt, &n, &n, &p, &tol, &rank,qraux, pivot, work);

    if (rank != p) singular=NOTFOUND;
    /* error("Singular matrix in qr_solve\n"); */

    coeft=dvector((*n1)*(*n1));

    /** Copy the matrix by column **/
    yt = dvector((*n1)*(*n1));
    for(i=0;i<*n1;i++)
      for(j=0;j<*n1;j++)
	yt[i*(*n1)+j]=y[j][i];

    F77_CALL(dqrcf)(xt, &n, &rank, qraux,yt, &n, coeft, &info);

    /** Put back into a matrix **/
    for(i=0;i<*n1;i++)
      for(j=0;j<*n1;j++)
	coef[j][i]=coeft[i*(*n1)+j];
    
    free(qraux);
    free(pivot);
    free(work);	
    free(xt);
    free(yt);
    free(coeft);
}


/* Thanks to Raphael Gottardo */
int inverse(double **mat1, int *n ,double **res)
{

  /** QR decomposition solve mat1*x=I **/
  /** res contains the result **/
  int i, flag=FOUND, nn=*n;
  double **iden;

  iden=dmatrix(nn,nn);    /* iden=dmatrix(nn,nn);   */
  for(i=0;i<*n;i++)
    iden[i][i]=1.0;
  qr_solve(mat1, n, iden, res, flag);
  free_dmatrix(iden,*n);
  return(flag);
}

/* vectors is non zero if want eigen vectors returned.
   length(EValues) = n
   dim(EVectors) = n,n
*/
int sym_eigen(double **A, int n, int vectorsflag, double *EValues, double **EVectors)
{
   int err=0,i=0,j=0;
   double *vA = dvector(n*n);
   double *vEVectors = dvector(n*n);
   double *l=dvector(n);
   double *fv1 = dvector(n);
   double *fv2 = dvector(n);
   double *wr = dvector(n);
   double *wi= dvector(n);
   double *z = dvector(n*n);
   int *iv1 = ivector(n);

   /** Copy the matrix by column **/
   /* make A into a vector to pass in*/
   for(j=0;j<n;j++){
     for(i=0;i<n;i++){
       vA[(n*j)+i] = A[i][j];
     }
   }

 
   /*  int F77_NAME(rs)(int *nm, int *n, double *a, double *w,
       int *matz, double *z, double *fv1, double *fv2, int *ierr) */
   F77_NAME(rs)(&n, &n, vA, l, &vectorsflag, vEVectors, fv1, fv2, &err);

   for(i=0;i<n;i++){
     EValues[i] = l[n-1-i];
   }  
   
   /* put A and EVectors into matrices to return */
   for(i=0;i<n;i++){
     for(j=0;j<n;j++){
       EVectors[i][j] = vEVectors[i + (n-(j+1))*n];
     }
   } 

   free_dvector(vA);
   free_dvector(vEVectors);
   free_dvector(l);
   free_dvector(fv1);
   free_dvector(fv2);
   free_dvector(wr);
   free_dvector(wi);
   free_dvector(z);
   free_ivector(iv1);
   return 0;
}



