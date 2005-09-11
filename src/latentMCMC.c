/****************************************************************************/
/*  Author: Susan Shortreed, susanms@stat.washington.edu                    */
/*  Purpose: main functions for parameter estimation for the latent         */
/*           space model proposed by  Peter D. Hoff, Adrian E. Raftery      */
/*           and Mark S. Handcock in                                        */
/*           "Latent Space Approaches to Social Network Analysis"           */
/*           All of this code is for an R function to be incorporated       */
/*           into the R ERGMM package.                                       */
/****************************************************************************/
/*####                 General Notes good to know                      #####*/
/*#  Z is the matrix of positions                                           */
/*#  k is the  dimension of the latent space                                */
/*#  g is the number of actors/nodes in the graph                           */
/*#                                                                         */
/****************************************************************************/
#include <math.h>
#include "R.h"
#include "Rmath.h"
/*#include "latentfunctionsb.c"*/
#include "latentUtil.h"
#include "latentfunctionsb.h"

#define TOLERANCE 1e-10
#define MAXIT 10000
#define TMAX 10
#define TRACE -1
#define TEMP 10
#define OLD 0
#define NEW 1
#define FOUND 0
#define NOTFOUND 1

void ergmm_latent(int* heads, int* tails, 
		 int n_edges, int g,
		 int MCMCSampleSize, int burnin, 
		 int interval, int k, 
		 double z_prior_mu, double z_prior_sd,
		 double b_prior_mu, double b_prior_sd, 
		 double **Z_mle, double zdelta,
		 double *vZ_post, double *Z_rate,
		 double *beta_mle, double bdelta,
		 double *Beta, double *B_rate, 
		 int p, int dir, double ***X,
		 double *Llike);

void MCMC_latent_wrapper( int *heads, int *tails, 
			  int *n_edges,	int *n,
			  int *MCMCSampleSize, int *burnin, 
			  int *interval, int *dimSpace, 
			  double *z_prior_mu, double *z_prior_sd,
			  double *b_prior_mu, double *b_prior_sd,
			  double* vZ_mle, double *Z_delta,
			  double *vZ_post, double *Z_rate,
			  double *beta_mle, double *beta_delta,
			  double *Beta, double *B_rate, 
			  int *p, int *dir, double *vX,
			  double *Llike
);


/* Here is the wrapper function we will call from the ergmm R code */
/* NOTE THE RANDOM SEED MUST BE SET IN R*/
void MCMC_latent_wrapper( int *heads, int *tails, int *n_edges,	int *n,
			  int *MCMCSampleSize, int *burnin, 
			  int *interval, int *dimSpace, 
			  double *z_prior_mu, double *z_prior_sd,
			  double *b_prior_mu, double *b_prior_sd,
			  double* vZ_mle, double *zdelta,
			  double *vZ_post, double *Z_rate,
			  double *beta_mle, double *bdelta,
			  double *Beta, double *B_rate, 
			  int *p, int *dir, double *vX,
			  double *Llike
)
{
  int i=0,j=0,k=0;
  /*   int nsample = (((*MCMCSampleSize)-(*burnin))/(*interval)); */
  double **Z_mle = dmatrix(*n,*dimSpace);
  double ***X = malloc( (*p) * sizeof(double**) );

  /*  set up all of the covariate matrices if covariates are involed  */
  /*  if p=0 (ie no covariates then these next 2 loops will do nothing) */
  for(i=0;i<*p;i++){
    X[i] = dmatrix((*n),(*n));
  } 
  for(k=0;k<*p;k++){
    for(i=0;i<*n;i++){
      for(j=0;j<*n;j++){
	X[k][i][j] = vX[ k*(*n)*(*n) + i*(*n) + j ];
      }
    }
  }

  /* Set Z up as a matrix, starting at the mle*/
  for(j=0;j<*dimSpace;j++){
    for(i=0;i<*n;i++){
      Z_mle[i][j] = vZ_mle[i+j*(*n)];
    }
  }

  /* R function enabling uniform RNG */
  GetRNGstate(); 
  ergmm_latent(heads, tails,*n_edges, *n, *MCMCSampleSize, *burnin, 
	      *interval, *dimSpace,
	      *z_prior_mu, *z_prior_sd,
	      *b_prior_mu, *b_prior_sd,
	      Z_mle, *zdelta, vZ_post, Z_rate, beta_mle, *bdelta,
	      Beta, B_rate, *p, *dir, X, Llike);
  PutRNGstate();
  
  /* free the memory used*/
  for(i=0;i<*p;i++){
    free_dmatrix(X[i],*n);
  }
  free(X);
  free_dmatrix(Z_mle,*n);
  return;
}
void ergmm_latent(int* heads, int* tails, 
		 int n_edges, int g,
		 int MCMCSampleSize, int burnin, 
		 int interval, int k, 
		 double z_prior_mu, double z_prior_sd,
		 double b_prior_mu, double b_prior_sd,
		 double **Z_mle, double zdelta,
		 double *vZ_post, double *Z_rate,
		 double *beta, double bdelta,
		 double *Beta, double *B_rate, 
		 int p, int dir, double ***X,
		 double *Llike)
{
   double **Z, **Znew, **pZ;
   /*  double *avZ, *param; */
   /* double *vZ, *svZ, *NRvZ, *AvZ; */
   double lik = 0.0;
   int n_accept_z=0, n_accept_b=0, top=0;
   int i=0,j=0,l=0, n_sample=0;
   /*    int **D; */
   int tempINT;

   n_sample = (int)((MCMCSampleSize-burnin)/interval);

   Z=dmatrix(g,k);
   copy_dmatrix(Z_mle,Z,g,k); 
   Znew = dmatrix(g,k);
   pZ=dmatrix(g,k);

   /*    printf("Here is the initial Z:\n"); */
   /*    print_dmatrix(Z,g,k,stdout); */
    
   lik = loglike_y(heads,tails,n_edges,g,Z,k,beta,p,dir,X);

   for(i=0;i<MCMCSampleSize;i++){

     /* update Z */
     /*      Rprintf("Testing if new.\n"); */
     tempINT = Z_up(heads,tails,n_edges,Z,zdelta,z_prior_mu,z_prior_sd,
		    g,k,Znew,beta,p,dir,&lik,X);
     /* Rprintf("lik = %f %d\n",lik,tempINT); */
     /*   Rprintf("We found Z_up = %d whereas NEW = %d\n",tempINT,NEW); */
     if(tempINT == NEW){
	 n_accept_z++;
	 copy_dmatrix(Znew,Z,g,k);	 
     }
     init_dmatrix(Znew,g,k,0.0);	     

     /* update beta given this new value of Z*/
     /* and conditioned on everything else*/
     /*      lik=loglike_y(heads,tails,n_edges,g,Z,k,beta,p,dir,X); */
     if( beta_up(heads,tails,n_edges,Z,g,k,b_prior_mu,
 	   b_prior_sd,&lik,beta,p,dir,X,bdelta) == NEW )
     {
      /*   Rprintf("lik = %f\n",lik); */
	   n_accept_b++;
     } 
     /* every interval after the burnin save the results */
     if( (i-burnin)>=0 && (i-burnin)%interval == 0){       
       top = ((i-burnin)/interval); 
       for(j=0;j<p;j++){
         Beta[j*n_sample+top] = beta[j];
       }
       /* we need to do procrustes on sampled Z */
/*      	 	 printf("Here is the untransformed Z:\n"); */
/*        	 	 print_dmatrix(Z,g,k,stdout); */
/*        if(procr_transform(Z,Z_mle,g,k,pZ) == FOUND){ */
/*          init_dmatrix(Z,g,k,0.0);  */
/* 	 copy_dmatrix(pZ,Z,g,k);	  */
/* 	 init_dmatrix(pZ,g,k,0.0);  */
/*        } */
       Llike[top] = lik; 
/*         Rprintf("lik = %f, Llike %f, top %d\n",lik,Llike[top],top); */
       for(j=0;j<k;j++){ 
	 for(l=0;l<g;l++){  
	   vZ_post[top*g*k+j*g+l] = Z[l][j];  
	 } 
       } 
       Z_rate[top] = (double) ((double)n_accept_z)/((double)interval); 
/*       printf("Z_rate %lf\n",Z_rate[top]); */
       B_rate[top] = (double) ((double)n_accept_b)/((double)interval);
       n_accept_z=0; 
       n_accept_b=0; 
     }    

    }

   /*Clean up memory mess*/ 
   free_dmatrix(Znew,g);
   free_dmatrix(Z,g);
   free_dmatrix(pZ,g);

   return;
  
}
