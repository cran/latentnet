/****************************************************************************/
/*  Author: Susan Shortreed, susanms@stat.washington.edu                    */
/*  Purpose: support functions for latentMCMC.c which uses the latent       */
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

#include <R.h>
#include <Rmath.h>
#include <math.h> 
#include <stdio.h>
#include <stdlib.h>
#include "latentUtil.h"
#include "latentfunctionsb.h"

#define OLD 0
#define NEW 1
#define FOUND 0
#define NOTFOUND 1


/* log probability of the graph, using following equation 
   logP(Y|n) = sum( eta[i,j]*Y[i,j] - log( 1+exp(eta[i,j]) ) )
   where eta = logodds(Y[i,j]=1|Z,a,b) = a - |Z[,i] - Z[,j]|
   Y is sociomatrix
   Z are latent positions */
double loglike_yplot(int *heads, int *tails, int n_edges, int g, double **Z, 
                 int k, double *beta, int dir, int nsubsamp)
{
  double llk=0.0;
  double etaij, tempijl, templ, u, sampfrac;
  int i, j, l,loop;

  /* actually calculating likelihood of y which is
   logP(Y|eta) = sum( eta[i][j]*y[i][j] - log(1+exp(eta[i][j]))  ) + g*log2
   the last term is due to the diag elements summed in R */

  for(i=0;i<n_edges;i++){
    tempijl=0.0;
    for(l=0;l<k;l++){
     templ = Z[heads[i]][l] - Z[tails[i]][l];
     tempijl += templ * templ;
    }
    etaij = beta[0] - sqrt(tempijl);
    llk += etaij;
  }

  sampfrac = (double)((g*(g-1))/2) / (double)nsubsamp;

  if(dir==NEW){
   if(nsubsamp >= (g*(g-1))/2){
    for(i=0;i<g;i++){
     for(j=0;j<i;j++){
       tempijl=0.0;
       for(l=0;l<k;l++){
        templ = Z[i][l] - Z[j][l];
        tempijl += templ * templ;
       }
       etaij = beta[0] - sqrt(tempijl);
       llk -= 2.0*log( 1.0+exp(etaij) );
     }
    }
   }else{
    for(loop=0;loop<nsubsamp;loop++)
    {
      u = runif(0.0,1.0)*(double)(g*(g-1));
      i = floor(u / (g-1) );
      j = floor(u - i*(g-1));
      if(j >= i){j++;}
/*     Rprintf("i %d j %d\n",i, j); */
       tempijl=0.0;
       for(l=0;l<k;l++){
        templ = Z[i][l] - Z[j][l];
        tempijl += templ * templ;
       }
       etaij = beta[0] - sqrt(tempijl);
       llk -= 2.0*log( 1.0+exp(etaij) );
    }
    llk -= log(2.0*sampfrac);
   }
  }else{
   if(nsubsamp >= (g*(g-1))/2 ){
    for(i=0;i<g;i++){
     for(j=0;j<i;j++){
       tempijl=0.0;
       for(l=0;l<k;l++){
        templ = Z[i][l] - Z[j][l];
        tempijl += templ * templ;
       }
       etaij = beta[0] - sqrt(tempijl);
       llk -= log( 1.0+exp(etaij) );
     }
    }
   }else{
    for(loop=0;loop<nsubsamp;loop++)
    {
      u = runif(0.0,1.0)*(double)(g*(g-1));
      i = floor(u / (g-1) );
      j = floor(u - i*(g-1));
      if(j >= i){j++;}
       tempijl=0.0;
       for(l=0;l<k;l++){
        templ = Z[i][l] - Z[j][l];
        tempijl += templ * templ;
       }
       etaij = beta[0] - sqrt(tempijl);
       llk -= log( 1.0+exp(etaij) );
    }
    llk -= log(sampfrac);
   }
  }

  return(llk);
}


/* update Z
   returns the updated value of the position matrix which is 
   is stored in Z.  ie it replaces the Z passed in. */
/*   
   head and tails form an edgelist
   n_edges is the number of edges in the graph (length of both heads and tails)
   Z are latent positions
   zdelta is the std. deviation on the proposal distr (Normal) for the Z[i,j]
   prior_mu and prior_std are the prior mean and standard deviation on the positions    
   g is the number of actors
   k is the number of dimensions
   Znew is the matrix for the propsed Z 
   beta is the covariate coefficients in the likelihood
   p is the number of covariates included in the model
   X is a vector of covariate difference matrices
*/
  int Z_upplot(int *heads, int *tails, int n_edges, double **Z, double zdelta, 
	     double prior_mu, double prior_sd, int g, int k,
	     double **Znew, double *beta, int dir, 
	     double *llk, int nsubsamp)
{
  double llknew=0.0, llkold, lr=0.0, temp;
  int i, j;

  /*   Rprintf("zdelta = %lf\n",zdelta); */

  /* perturb the old z values with a N(0,zdelta) distribution */
  for(i=0;i<g;i++){
    for(j=0;j<k;j++){
      Znew[i][j] = Z[i][j] + rnorm(0.0,zdelta);
    }
  }
  /*    printf("Here is the testing Z:\n"); */
  /*    print_dmatrix(Znew,g,k,stdout); */

  /* calculate the loglikelihoods */
  llknew = loglike_yplot(heads,tails,n_edges,g,Znew,k,beta,dir,nsubsamp);
  llkold = *llk;
  /*   printf("Old Likelihood = %lf and New Likelihood = %lf:\n",llkold,llknew); */
  /* calculate the acceptance ratio */
  lr = llknew-llkold;

  /* this is the "Q" part */
  for(i=0;i<g;i++){
    for(j=0;j<k;j++){
      lr += dnorm(Znew[i][j],prior_mu,prior_sd,1);
      lr -= dnorm(Z[i][j],   prior_mu,prior_sd,1);
    }
  }  

  /* accept with probability r */
  temp = runif(0.0,1.0); 
  /*   printf("runif(0,1) = %lf\n",temp); */
  if( temp < exp(lr) ){
    /*     printf("we accepted it!\n"); */
    *llk = llknew;
    return(NEW);
  }
  else{
    /*     printf("we rejected it!\n"); */
    *llk = llkold;
    return(OLD);
  }
}

/* updates beta
   returns the loglikelihood for the updated value of beta, the updated value of beta 
   is stored in beta.  ie it replaces the beta passed in. */
/*   
   head and tails form an edgelist
   n_edges is the number of edges in the graph (length of both heads and tails)
   Z are latent positions
   g is the number of actors
   k is the number of dimensions
   prior_mu and prior_sd are the parameters for the Normal prior on beta
   llk is for the value of the likelihood
   beta is the covariate coefficients in the likelihood
   p is the number of covariates included in the model
   X is a vector of covariate difference matrices
   bdelta is the standard deviaton of the proposal distribution which is a 
           multivariate normal centered around the previous beta value
*/
int beta_upplot(int *heads, int *tails, int n_edges, double **Z, int g, int k, 
            double prior_mu, double prior_sd, double *llk,
            double *beta, int dir, 
	    double bdelta,int nsubsamp)
{
  double llknew=0.0, llkold, lr=0.0;
  double *betanew = dvector(1);

  betanew[0] = beta[0] + rnorm(0.0,bdelta);

  llknew = loglike_yplot(heads,tails,n_edges,g,Z,k,betanew,dir,nsubsamp);
  llkold = *llk;
  /* calculate the acceptance ratio */
  lr = llknew-llkold;
  /* this is the "Q" part */
  lr += dnorm(betanew[0],prior_mu,prior_sd,1);
  lr -= dnorm(beta[0],   prior_mu,prior_sd,1);  
  
  /* accept with probability r */
  if( runif(0.0,1.0) < exp(lr) ){
    beta[0] = betanew[0];
    *llk = llknew; 
    free_dvector(betanew);
    return(NEW);
  }
  else{
    *llk = llkold;
    free_dvector(betanew);
    return(OLD);
  }
}
