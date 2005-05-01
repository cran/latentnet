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
double loglike_y(int *heads, int *tails, int n_edges, int g, double **Z, 
                 int k, double *beta, int p, int dir, double ***X)
{
  double llk=0.0;
  double etaij, tempijl, templ;
  int i, j, l;
  double **cov = dmatrix(g,g);
  double **gg_helper = dmatrix(g,g);

  /* get covariate matrix ready to add on, if no covariates nothing done */
  for(i=0;i<p;i++){
    dscaler_times_matrix(beta[i],X[i],g,g,gg_helper);
    dmatrix_addition(cov,g,g,gg_helper);
  }

  /* actually calculating likelihood of y which is
   logP(Y|eta) = sum( eta[i][j]*y[i][j] - log(1+exp(eta[i][j]))  ) + g*log2
   the last term is due to the diag elements summed in R */

  for(i=0;i<n_edges;i++){
    tempijl=0.0;
    for(l=0;l<k;l++){
     templ = Z[heads[i]][l] - Z[tails[i]][l];
     tempijl += templ * templ;
    }
    etaij = cov[heads[i]][tails[i]] - sqrt(tempijl);
    llk += etaij;
  }

  if(dir==NEW){
   for(i=0;i<g;i++){
    for(j=0;j<i;j++){
       tempijl=0.0;
       for(l=0;l<k;l++){
        templ = Z[i][l] - Z[j][l];
        tempijl += templ * templ;
       }
       etaij = cov[i][j] - sqrt(tempijl);
       llk -= log( 1.0+exp(etaij) );
       etaij = cov[j][i] - sqrt(tempijl);
       llk -= log( 1.0+exp(etaij) );
    }
   }
  }else{
   for(i=0;i<g;i++){
    for(j=0;j<i;j++){
       tempijl=0.0;
       for(l=0;l<k;l++){
        templ = Z[i][l] - Z[j][l];
        tempijl += templ * templ;
       }
       etaij = cov[i][j] - sqrt(tempijl);
       llk -= log( 1.0+exp(etaij) );
    }
   }
  }

  /* Clean the memory mess */
  free_dmatrix(cov,g);
  free_dmatrix(gg_helper,g);
  return(llk);
}

/*  Gives the negative distance between nodes
    for dim k: d[i,j] <- sqrt( (x1i-x2j)^2 + ..... + (xki-xkj)^2 ) */
/*
  Z is position matrix
  g is number of actors
  dim is number of dimensions of latent space
  dz is the distance matrix that is returned
*/
/*  dZ must be a g by g matrix */
void neg_dist_Z(double **Z, int g, int dim, double **dZ)
{
  double tempk, tempijk;
  int i,j,k;
  
  for(i=0;i<g;i++){
    for(j=i+1;j<g;j++){
      tempijk=0.0;
      for(k=0;k<dim;k++){
        tempk = Z[i][k] - Z[j][k];
        tempijk += tempk * tempk; 
      }
      dZ[i][j] = - sqrt(tempijk);
      dZ[j][i] = dZ[i][j];
    }
  }
  return;
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
int Z_up(int *heads, int *tails, int n_edges, double **Z, double zdelta, 
         double prior_mu, double prior_sd, int g, int k,
         double **Znew, double *beta, int p, int dir, 
	 double *llk, double ***X)
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
  llknew = loglike_y(heads,tails,n_edges,g,Znew,k,beta,p,dir,X);
  llkold = *llk;
  /*  llkold = loglike_y(heads,tails,n_edges,g,Z,   k,beta,p,dir,X); */
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

  /*   printf("prior_mu = %lf and prior_sd = %lf so that we get lr = %lf\n",prior_mu,prior_sd,lr); */

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
int beta_up(int *heads, int *tails, int n_edges, double **Z, int g, int k, 
            double prior_mu, double prior_sd, double *llk,
            double *beta, int p, int dir, double ***X,
	    double bdelta)
{
  int i=0;
  double llknew=0.0, llkold, lr=0.0;
  double *betanew = dvector(p);

  for(i=0;i<p;i++){
    betanew[i] = beta[i] + rnorm(0.0,bdelta);
  }

  llknew = loglike_y(heads,tails,n_edges,g,Z,k,betanew,p,dir,X);
  llkold = *llk;
  /*  llkold = loglike_y(heads,tails,n_edges,g,Z,k,beta,   p,dir,X); */
  /* calculate the acceptance ratio */
  lr = llknew-llkold;
  /* this is the "Q" part */
  for(i=0;i<p;i++){
    lr += dnorm(betanew[i],prior_mu,prior_sd,1);
    lr -= dnorm(beta[i],   prior_mu,prior_sd,1);  
  }
  
  /* accept with probability r */
  if( runif(0.0,1.0) < exp(lr) ){
    for(i=0;i<p;i++){
      beta[i] = betanew[i];
    }
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


/*
 Finds Procrustes transform (denoted Z* below).  
 Gives the translation, rotation, reflection of Z closet to Zo
 The formula is Z* = Z[(Z'ZoZo'Z)^(-1/2)]Z'Zo
     when Z and Zo are centered around the same origin
     and Z, Zo, Z* are n x k martices (transpose of the paper)  
 k is the dimension of the space that we are working in.
 note this is slightly different then the paper,
 in the paper Z is k by n, while her is n by k
*/
int procr_transform(double **Z, double **Zo, int g, int k, double **pZ)
{
  
  int i=0,j=0;
  double **A, **tZ, **tZo, ** Ahalf, **AhalfInv, **tptrans;
  double **eAvectors, **eADvalues, **teAvectors;
  double *avZ, *avZo, *eAvalues;
  double **kk_helper, **gg_helper, **gk_helper, **kg_helper, **kk2_helper;
  /*   double *tempZ, *tempZo, *temp ; */
  kg_helper=dmatrix(k,g);
  gk_helper=dmatrix(g,k);
  kk_helper=dmatrix(k,k);
  gg_helper=dmatrix(g,g);
  kk2_helper=dmatrix(k,k);

  /*  Center both Z and Zo around the mean of Zo    */
  /*  First centers Z around the origin  */
  avZ = dvector(k);
  avZo = dvector(k);
  
  for(j=0;j<k;j++){
    for(i=0;i<g;i++){
      avZ[j] += Z[i][j] / g;
      avZo[j] += Zo[i][j] / g;
    }
  }

  /*   subtract the averages */ 
  for(j=0;j<k;j++){
    for(i=0;i<g;i++){
      Z[i][j] -= (avZ[j] - avZo[j]);
    }
  }

  /* Compute A = tZ*Zo*tZo*Z*/
  A=dmatrix(k,k);
  tZ=dmatrix(k,g);
  tZo=dmatrix(k,g);
  t(Z,g,k,tZ);
  t(Zo,g,k,tZo);

  dmatrix_multiply(tZ,k,g,Zo,k,kk_helper);
  dmatrix_multiply(kk_helper,k,k,tZo,g,kg_helper);
  dmatrix_multiply(kg_helper,k,g,Z,k,A);
  init_dmatrix(kk_helper,k,k,0.0);
  init_dmatrix(kg_helper,k,g,0.0);

  /* Compute sqrt(A) */
  eAvectors=dmatrix(k,k);
  eAvalues=dvector(k);
  sym_eigen(A,k,1,eAvalues,eAvectors);
  eADvalues=dmatrix(k,k);
  for(i=0;i<k;i++){
    eADvalues[i][i] = sqrt(eAvalues[i]);
  }
  teAvectors=dmatrix(k,k);
  t(eAvectors,k,k,teAvectors);  

  Ahalf=dmatrix(k,k);
  dmatrix_multiply(eAvectors,k,k,eADvalues,k,kk_helper);
  dmatrix_multiply(kk_helper,k,k,teAvectors,k,Ahalf);
  init_dmatrix(kk_helper,k,k,0.0);

  /* Now compute the inverse */
  AhalfInv=dmatrix(k,k);
  if(inverse(Ahalf,&k,AhalfInv) != FOUND) return NOTFOUND;

  /*Now compute t(Zo)*Z*AhalfInv*tZ*/
  tptrans=dmatrix(k,g);
  dmatrix_multiply(tZo,k,g,Z,k,kk_helper);
  dmatrix_multiply(kk_helper,k,k,AhalfInv,k,kk2_helper);
  dmatrix_multiply(kk2_helper,k,k,tZ,g,tptrans);
  t(tptrans,k,g,pZ);

  free_dmatrix(A,k);
  free_dmatrix(tZ,k);
  free_dmatrix(tZo,k);
  free_dmatrix(gk_helper,g);
  free_dmatrix(kk_helper,k);
  free_dmatrix(gg_helper,g);
  free_dmatrix(kk2_helper,k);
  free_dmatrix(kg_helper,k);
  free_dmatrix(eAvectors,k);
  free_dmatrix(eADvalues,k);
  free_dmatrix(teAvectors,k);
  free_dmatrix(Ahalf,k);
  free_dmatrix(AhalfInv,k);
  free_dmatrix(tptrans,k);
  free_dvector(avZ);
  free_dvector(avZo);
  free_dvector(eAvalues);
    
  return(FOUND);
}
