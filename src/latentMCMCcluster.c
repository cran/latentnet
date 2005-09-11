/****************************************************************************/
/*  Original Author: Susan Shortreed, susanms@stat.washington.edu           */
/*  Updated by: Jeremy Tantrum, tantrum@stat.washington.edu                 */
/*  Purpose: main functions for parameter estimation model 2                */
/*           All of this code is for an R function which is incorporated    */
/*           into the R ERGMMM package.                                       */
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
#include "latentUtil.h"
#include "latentfunctionsbcluster.h"
#include "latentfunctionsb.h"
#include "wishart.h"
#include "mvnorm.h"

#define TOLERANCE 1e-10
#define MAXIT 10000
#define TMAX 10
#define TRACE -1
#define TEMP 10
#define OLD 0
#define NEW 1
#define FOUND 0
#define NOTFOUND 1



void ergmm_latent2(int* heads, int* tails, 
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
		  double *Llike, int ng,
		  double *epsilon, double *mu, double *Sigma, int *Ki,
		  long *KiList, double *muList, double *SigmaList,
		  double Sigprior, double muSigprior, double dirprior, 
		  double alphaprior, int chisqprop, double thetaprop);

void MCMC_latent2_wrapper( int *heads, int *tails, 
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
			   double *Llike, int *ng,
			   double *epsilon,double *mu,double *Sigma,int *Ki,
			   long *KiList, double *muList, double *SigmaList,
			   double *Sigprior, double *muSigprior, 
			   double *dirprior, double *alphaprior,
			   int *chisqprop, double *thetaprop);


/* Here is the wrapper function we will call from the ergmm R code */
/* NOTE THE RANDOM SEED MUST BE SET IN R*/
void MCMC_latent2_wrapper( int *heads, int *tails, int *n_edges, int *n,
			   int *MCMCSampleSize, int *burnin, 
			   int *interval, int *dimSpace, 
			   double *z_prior_mu, double *z_prior_sd,
			   double *b_prior_mu, double *b_prior_sd,
			   double* vZ_mle, double *zdelta,
			   double *vZ_post, double *Z_rate,
			   double *beta_mle, double *bdelta,
			   double *Beta, double *B_rate, 
			   int *p, int *dir, double *vX,
			   double *Llike, int *ng,
			   double *epsilon, double *mu,double *Sigma,int *Ki,
			   long *KiList, double *muList, double *SigmaList,
			   double *Sigprior, double *muSigprior, 
			   double *dirprior, double *alphaprior,
			   int *chisqprop, double *thetaprop)
{
  int i=0,j=0,k=0;
  double **Z_mle = dmatrix(*n,*dimSpace);
  double ***X = malloc( (*p) * sizeof(double**) );

  /*  set up all of the covariate matrices if covariates are involed  */
  /*  if p=0 (ie no covariates then these next 2 loops will do nothing) */
  /*  */
  /*  Rprintf("Entered into MCMC_latent2_wrapper\n"); */

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

  /*   for(i=0;i<(*dimSpace);i++) */
  /*     for(j=0;j<(*ng);j++) */
  /*     { */
  /*       mumat[j][i] = mu[ i* (*ng) + j ]; */
  /*     } */

  /* Set Z up as a matrix, starting at the mle*/
  for(j=0;j<*dimSpace;j++){
    for(i=0;i<*n;i++){
      Z_mle[i][j] = vZ_mle[i+j*(*n)];
    }
  }

  /* R function enabling uniform RNG */
  GetRNGstate();
 
  ergmm_latent2(heads, tails,*n_edges, *n, *MCMCSampleSize, *burnin, 
	       *interval, *dimSpace,
	       *z_prior_mu, *z_prior_sd,
	       *b_prior_mu, *b_prior_sd,
	       Z_mle, *zdelta, vZ_post, Z_rate, beta_mle, *bdelta,
	       Beta, B_rate, *p, *dir, X, Llike,*ng,epsilon,mu,Sigma,Ki,
	       KiList, muList, SigmaList, *Sigprior, *muSigprior, *dirprior,
	       *alphaprior, *chisqprop, *thetaprop);
  PutRNGstate();
  
  /* free the memory used*/
  for(i=0;i<*p;i++){
    free_dmatrix(X[i],*n);
  }
  free(X);
  free_dmatrix(Z_mle,*n);
  return;
}

void ergmm_latent2(int* heads, int* tails, 
		  int n_edges, int g,
		  int MCMCSampleSize, int burnin, 
		  int interval, int dim, 
		  double z_prior_mu, double z_prior_sd,
		  double b_prior_mu, double b_prior_sd,
		  double **Z_mle, double zdelta,
		  double *vZ_post, double *Z_rate,
		  double *beta, double bdelta,
		  double *Beta, double *B_rate, 
		  int p, int dir, double ***X,
		  double *Llike, int ng,
		  double *epsilon, double *mu, double *Sigma, int *Ki,
		  long *KiList, double *muList, double *SigmaList, 
		  double Sigprior, double muSigprior, double dirprior,
		  double alphaprior, int chisqprop, double thetaprop)
{
  double **Z, **Znew, **pZ, *pK, temp, **mubar, Shat;
  double useSig, *vec1, *choldecomp, sum, Sbar, rms;
  double lik = 0.0;
  int n_accept_z=0, n_accept_b=0, top=0, i=0,j=0,l=0, n_sample=0, *n, tempINT;
  int loop1, mcmcloop;
  long longdim;

  /* memory for procr_transform2 */
  /* procrustes vars */
  double **A, **tZ, **tZo, ** Ahalf, **AhalfInv, **tptrans;
  double **eAvectors, **eADvalues, **teAvectors;
  double *avZ, *avZo, *eAvalues;
  double **kk_helper, **gg_helper, **gk_helper, **kg_helper, **kk2_helper;

  kg_helper=dmatrix(dim,g);
  gk_helper=dmatrix(g,dim);
  kk_helper=dmatrix(dim,dim);
  gg_helper=dmatrix(g,g);
  kk2_helper=dmatrix(dim,dim);
  avZ = dvector(dim);
  avZo = dvector(dim);
  A=dmatrix(dim,dim);
  tZ=dmatrix(dim,g);
  tZo=dmatrix(dim,g);
  eAvectors=dmatrix(dim,dim);
  eAvalues=dvector(dim);
  eADvalues=dmatrix(dim,dim);
  teAvectors=dmatrix(dim,dim);
  Ahalf=dmatrix(dim,dim);
  AhalfInv=dmatrix(dim,dim);
  tptrans=dmatrix(dim,g);


  n_sample = (int)((MCMCSampleSize-burnin)/interval);

  Z=dmatrix(g,dim);
  copy_dmatrix(Z_mle,Z,g,dim); 
  Znew = dmatrix(g,dim);
  pZ=dmatrix(g,dim);
  n = ivector(ng); /*  n[i] is the number of obs in group i - amoung ng groups */
  pK = dvector(ng);
  mubar = dmatrix(ng,dim);
  vec1 = dvector(dim);
  choldecomp = dvector(dim * (dim-1));


  for(i=0;i<g;i++)
    n[Ki[i] - 1]++;

  for(mcmcloop=0;mcmcloop<MCMCSampleSize;mcmcloop++){
    /*     Rprintf("iteration %3d - ",mcmcloop); */
    /* update Z */
    lik = loglike2_y(heads,tails,n_edges,g,Z,dim,beta,p,dir,X, mu, Sigma, Ki, ng);
    tempINT = Z_up2(heads,tails,n_edges,Z,zdelta,z_prior_mu,z_prior_sd,
		    g,dim,Znew,beta,p,dir,X,mu,Sigma,Ki,ng,&lik,
		    chisqprop, thetaprop, Z_mle, A, tZ, tZo, Ahalf, AhalfInv, 
		    tptrans, eAvectors, eADvalues, teAvectors, avZ, 
		    avZo, eAvalues, kk_helper, gg_helper, gk_helper, 
		    kg_helper, kk2_helper);

    if(tempINT == NEW){
      n_accept_z++;
      copy_dmatrix(Znew,Z,g,dim);
    }
    /*   init_dmatrix(Znew,g,dim,0.0); */

/*     if(procr_transform(Z,Z_mle,g,dim,pZ) == FOUND){ */
/*       init_dmatrix(Z,g,k,0.0);  */
/*       copy_dmatrix(pZ,Z,g,k);	  */
/*       init_dmatrix(pZ,g,k,0.0);  */
/*     } */

    rms = 1.0; /*  +++++++++++ REMOVE */


    /* update beta given this new value of Z*/
    /* and conditioned on everything else*/

    lik=loglike2_y(heads,tails,n_edges,g,Z,dim,beta,p,dir,X, mu, Sigma, Ki, ng);
    if( beta_up2(heads,tails,n_edges,Z,g,dim,b_prior_mu,
		 b_prior_sd,&lik,beta,p,dir,X,bdelta,
		 mu,Sigma,Ki,ng) == NEW )
    {
      n_accept_b++;
    } 

    longdim = (long)dim;

    /*-- Update the other parameters --*/
    /* Ki ~ p(Ki[i]==j).... */
/*     if( (mcmcloop-burnin)>=0)  */
/*       for(j=0;j<ng;j++) */
/* 	Rprintf("mu_%d = %1.4lf  %1.4lf   -- S = %1.4lf\n",j,mu[j],mu[j+ng],Sigma[j]); */

/*     if( (mcmcloop-burnin)>=0) */
/*     { */
/*       Rprintf("mu_%d = (%1.4lf, %1.4lf) (%1.4lf, %1.4lf) (%1.4lf, %1.4lf) (%1.4lf, %1.4lf)\n",mcmcloop,mu[0],mu[0+ng],mu[1],mu[1+ng],mu[2],mu[2+ng],mu[3],mu[3+ng]); */
/*       Rprintf("sigma_%d = (%1.4lf %1.4lf %1.4lf %1.4lf)\n",mcmcloop,Sigma[0],Sigma[1],Sigma[2],Sigma[3]); */
/*     } */
    for(i=0;i<g;i++)
    {
      sum = 0;
      for(j=0;j<ng;j++)
      {
/* 	if( (mcmcloop-burnin)>=0)  */
/* 	{ */
/* 	  Sbar = (Z[i][0]-mu[j])*(Z[i][0]-mu[j])/ (-2*Sigma[j]); */
/* 	  Sbar += (Z[i][1]-mu[j+ng])*(Z[i][1]-mu[j+ng])/ (-2*Sigma[j]); */
/* 	  Rprintf("x-mu=%1.4f - x-mu/-s=%1.4lf lik=%1.4f exp = %1.4lf eps = %1.4f\n",Z[i][0]-mu[j],Sbar,-0.5 * longdim * log(2 * M_PI * Sigma[j]) + Sbar,exp(-0.5 * longdim * log(2 * M_PI * Sigma[j]) + Sbar), epsilon[j]); */
/* 	} */

	sdlnorm(&longdim,ng,j,mu,Sigma,Z[i], &temp);
	if(j>0)
	  pK[j] = pK[j-1] + epsilon[j] * exp(temp);
	else 
	  pK[j] = epsilon[j] * exp(temp);
      }

      temp = runif(0.0,1.0);
      j = 0;
      while(pK[j]/pK[ng-1] < temp)
	j++;
/*       if( (mcmcloop-burnin)>=0) */
/*       { */
/* 	Rprintf("pk[%d]= ",i); */
/* 	for(loop1=0;loop1<ng;loop1++) */
/* 	  Rprintf("%1.4f ",pK[loop1]/pK[ng-1]); */
/* 	Rprintf(" ==> %d\n",j+1); */
/*       } */

      Ki[i] = j + 1;
    }

    for(i=0;i<ng;i++)
      n[i] = 0;

    for(i=0;i<g;i++)
      n[Ki[i] - 1]++;

    /* rdirichlet for epsilon */
    for(i=0;i<ng;i++){
      epsilon[i] = (double)(n[i]+dirprior);
    }
    rdirich(ng,epsilon);

    init_dmatrix(mubar,ng,dim,0.0);
    for(i=0;i<g;i++)
      for(j=0;j<dim;j++)
	mubar[Ki[i]-1][j] += Z[i][j]/n[Ki[i]-1];


    /*  invchisq for sigma */
    Sbar = 0.0;
    for(i=0;i<ng;i++)
    {
      Shat = 0.0;
      for(j=0;j<g;j++)
	if((Ki[j] - 1) == i)
	    for(loop1=0;loop1<dim;loop1++)
	      Shat += (Z[j][loop1] - mu[i + loop1 * ng]) * (Z[j][loop1] - mu[i + loop1 * ng]);
      /*             Rprintf("n = %d, s^2 = %1.4f   <-> s^2/n = %1.4f\n",n[i],Shat/dim,Shat/(dim*n[i])); */
      if(n[i]>2)
      {
	temp = rchisq(n[i] + alphaprior);
/* 	Sigma[i] = (alphaprior *Sigprior + (n[i]-1) *Shat/dim +  */
/* 		    muSigprior*n[i]/(muSigprior+n[i]) *  */
/* 		    (mubar[i][0]*mubar[i][0] + mubar[i][1]*mubar[i][1])/2) /  */
/* 	  (n[i] + alphaprior) / temp; */
        	/*  2D model only. */

	Sigma[i] = (Sigprior +  Shat/dim) /temp;
	/* Sigma[i] = (alphaprior *Sigprior +  n[i]*Shat/dim)/(alphaprior+n[i]) /temp; */


	/* 	Rprintf("%1.4f,  %1f,\n",Sigprior+Shat/2,n[i]+alphaprior); */
	/* 	Rprintf("Sigma[%d] = %1.4f   rchisq() = %1.4f\n",i,Shat/dim,temp);  */
	/* 	Rprintf("Sigma[%d] = %1.4f, rms*Sigprior = %1.4f, epsilon = %1.4f, rchisq() = %1.4f\n",i,Sigma[i],rms*Sigprior,n[i]+alphaprior,temp);  */

	Sbar += Sigma[i];
      }
    }
    for(i=0;i<ng;i++)
    {
      if(n[i] < 3)
	Sigma[i] = Sbar/(ng - 1);
    }

    /* mvrnorm for mumat */

    for(i=0;i<ng;i++)
    {
      for(j=0;j<dim;j++)
	mubar[i][j] = n[i] * mubar[i][j]/(n[i]+ 1/muSigprior);
      /* 	mubar[i][j] = n[i] * mubar[i][j]/(n[i]+ muSigprior); */
      /* mubar[i][j] = (n[i] * mubar[i][j]/Sigma[i])/(1/muSigprior + n[i]/Sigma[i]); */
      /*       useSig = Sigma[i]/(n[i] + muSigprior); */
      useSig = Sigma[i]/(n[i] + Sigma[i]/muSigprior);
      /*       useSig = 1/(n[i] / Sigma[i] + 1/muSigprior); */
      /*       Rprintf("mu[%d] ~ N( (%1.4f %1.4f), %1.4f)\n",i,mubar[i][0],mubar[i][1],useSig); */
      for(j=0;j<dim;j++)
	mu[i + j * ng] = rnorm(mubar[i][j],sqrt(useSig));
      /*       Rprintf("mu[%d][0] = %1.4f\n",i,mu[i]); */
    }

    /* every interval after the burnin save the results */
    if( (mcmcloop-burnin)>=0 && (mcmcloop-burnin)%interval == 0){       
      top = ((mcmcloop-burnin)/interval); 
      for(j=0;j<=p;j++){/* beta[p] */
	Beta[j*n_sample+top] = beta[j];
      }
      Llike[top] = lik; 
      for(j=0;j<dim;j++){ 
	for(l=0;l<g;l++){  
	  vZ_post[top*g*dim+j*g+l] = Z[l][j];  
	} 
      } 
      Z_rate[top] = (double) ((double)n_accept_z)/((double)interval); 
      B_rate[top] = (double) ((double)n_accept_b)/((double)interval);
      for(j=0;j<g;j++)
      {
	KiList[top + j * n_sample] = Ki[j];
      }
      for(j=0;j<ng;j++)
      {
	for(l=0;l<dim;l++)
	  muList[top + n_sample * (l + j*dim)] = mu[j + l * ng];
	SigmaList[top + n_sample * j] = Sigma[j];
      }
      n_accept_z=0; 
      n_accept_b=0; 
    }    

  }

  /*Clean up memory mess*/ 
  free_dvector(choldecomp);
  free_dvector(vec1);
  free_ivector(n);
  free_dvector(pK);
  free_dmatrix(mubar,ng);

  free_dmatrix(Znew,g);
  free_dmatrix(Z,g);
  free_dmatrix(pZ,g);

  /*  FREE procrustes stuff! */
  free_dmatrix(A,dim);
  free_dmatrix(tZ,dim);
  free_dmatrix(tZo,dim);
  free_dmatrix(gk_helper,g);
  free_dmatrix(kk_helper,dim);
  free_dmatrix(gg_helper,g);
  free_dmatrix(kk2_helper,dim);
  free_dmatrix(kg_helper,dim);
  free_dmatrix(eAvectors,dim);
  free_dmatrix(eADvalues,dim);
  free_dmatrix(teAvectors,dim);
  free_dmatrix(Ahalf,dim);
  free_dmatrix(AhalfInv,dim);
  free_dmatrix(tptrans,dim);
  free_dvector(avZ);
  free_dvector(avZo);
  free_dvector(eAvalues);
  /* end free procrustes stuff! */


  return;
  
}
