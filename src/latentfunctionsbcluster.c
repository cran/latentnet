/****************************************************************************/
/*  Original Author: Susan Shortreed, susanms@stat.washington.edu           */
/*  Updated by: Jeremy Tantrum, tantrum@stat.washington.edu                 */
/*  Purpose: support functions for parameter estimation model 2             */
/*           All of this code is for an R function which is incorporated    */
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
/* dhunterchange:  Need following #include statements */
#include "latentUtil.h"
#include "latentfunctionsbcluster.h"
#include "mvnorm.h"

#define OLD 0
#define NEW 1
#define FOUND 0
#define NOTFOUND 1


/* log probability of the graph, using following equation 
   logP(Y|n) = sum( eta[i,j]*Y[i,j] - log( 1+exp(eta[i,j]) ) )
   where eta = logodds(Y[i,j]=1|Z,a,b) = a - |Z[,i] - Z[,j]|
   Y is sociomatrix
   Z are latent positions */
double loglike2_y(int *heads, int *tails, int n_edges, int g, double **Z, 
		  int k, double *beta, int p, int dir, double ***X, double *mu,
		  double *Sigma, int *Ki, int ng)
{
  double llk=0.0;
  double etaij, temp, tempijl, templ;
  int i, j, l;
  long longdim = k;
  double **cov = dmatrix(g,g); 
  double **gg_helper = dmatrix(g,g); 

  /* get covariate matrix ready to add on, if no covariates nothing done */
  init_dmatrix(cov,g,g,0.0);
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
    /* MISS::if(!(heads[i]==missing1 && tails[i]==missing2)){ */
    /*     etaij = cov[heads[i]][tails[i]]*(1 - sqrt(tempijl)); */
    etaij = cov[heads[i]][tails[i]] - beta[p] * sqrt(tempijl);
    /*     if((heads[i]==5) || tails[i]==5) */
    /*       Rprintf("beta[p] = %1.4f   dij = %1.4f\n",beta[p],sqrt(tempijl)); */
    llk += etaij;
  }
  /*   Rprintf("llk = %1.4f   ",llk); */

  if(dir==NEW){
    for(i=0;i<g;i++){
      for(j=0;j<i;j++){
	if(i!=j){
	  tempijl=0.0;
	  for(l=0;l<k;l++){
	    templ = Z[i][l] - Z[j][l];
	    tempijl += templ * templ;
	  }
	  /* MISS::if(!(heads[i]==missing1 && tails[i]==missing2)){ */
	  /* 	etaij = cov[i][j]*(1 - sqrt(tempijl)); */
	  etaij = cov[i][j] - beta[p] * sqrt(tempijl);
	  llk -= log( 1.0+exp(etaij) );
	  etaij = cov[j][i] - beta[p] * sqrt(tempijl);
	  llk -= log( 1.0+exp(etaij) );
	}
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
	etaij = cov[i][j] - beta[p] * sqrt(tempijl);
	llk -= log( 1.0+exp(etaij) );
      }
    }
  }

  /*   Rprintf("   llk = %1.4f\n",llk); */

  /*add in part of likelihood due to z_i's being random and not
   *parameters */
  for(i=0;i<g;i++){
    sdlnorm(&longdim,ng,Ki[i]-1,mu,Sigma,Z[i], &temp);
    llk += temp;
    /*     if(i==13) */
    /*       Rprintf("dlnorm = %1.4g  mu = (%1.4f,%1.4f) Sigma=%1.4f,  Z = (%1.4f, %1.4f)\n",temp,mu[Ki[i]-1],mu[Ki[i]-1+ng],Sigma[Ki[i]-1], Z[i][0], Z[i][1]); */
  }
  /*   Rprintf("temp = %1.4g   and  llk = %1.4g\n",temp,llk); */

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
void neg_dist_Z2(double **Z, int g, int dim, double **dZ)
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
int Z_up2(int *heads, int *tails, int n_edges, double **Z, double zdelta, 
	  double prior_mu, double prior_sd, int g, int k,
	  double **Znew, double *beta, int p, int dir, double ***X, 
	  double *mu,double *Sigma, int *Ki, int ng, double *llk, 
	  int chisqprop, double thetaprop, double **Z_mle,
	  double **A, double **tZ, double **tZo, double ** Ahalf, 
	  double **AhalfInv, double **tptrans, double **eAvectors, 
	  double **eADvalues, double **teAvectors, double *avZ, 
	  double *avZo, double *eAvalues, double **kk_helper, 
	  double **gg_helper, double **gk_helper, 
	  double **kg_helper, double **kk2_helper)
{
  /*   double llknew=0.0, llkold=0.0, lr=0.0, temp, *rn; */
  double llknew=0.0, llkold=0.0, lr=0.0, temp, usesigma;
  double **pZ, rms;
  int i, j,loop1, change=0;

  /* procrustes vars */
/*   double **A, **tZ, **tZo, ** Ahalf, **AhalfInv, **tptrans; */
/*   double **eAvectors, **eADvalues, **teAvectors; */
/*   double *avZ, *avZo, *eAvalues; */
/*   double **kk_helper, **gg_helper, **gk_helper, **kg_helper, **kk2_helper; */

/*   kg_helper=dmatrix(k,g); */
/*   gk_helper=dmatrix(g,k); */
/*   kk_helper=dmatrix(k,k); */
/*   gg_helper=dmatrix(g,g); */
/*   kk2_helper=dmatrix(k,k); */
/*   avZ = dvector(k); */
/*   avZo = dvector(k); */
/*   A=dmatrix(k,k); */
/*   tZ=dmatrix(k,g); */
/*   tZo=dmatrix(k,g); */
/*   eAvectors=dmatrix(k,k); */
/*   eAvalues=dvector(k); */
/*   eADvalues=dmatrix(k,k); */
/*   teAvectors=dmatrix(k,k); */
/*   Ahalf=dmatrix(k,k); */
/*   AhalfInv=dmatrix(k,k); */
/*   tptrans=dmatrix(k,g); */


  pZ=dmatrix(g,k);

  /* perturb the old z values with a N(mu_j,Sigma_j) distribution */
  /* g is the number of nodes */

  copy_dmatrix(Z,Znew,g,k);
  /*   Rprintf("About to start finding new Z\n"); */
  for(i=0;i<g;i++){
/*     rsq = Z[i][0] * Z[i][0] + Z[i][1] * Z[i][1]; */
/*     theta = asin(Z[i][1]/sqrt(rsq)); */
/*     if(Z[i][0] < 0) theta = M_PI - theta; */
/*     thetanew = rnorm(theta, thetaprop); */
/*     rsqnew = rchisq(chisqprop); */
/*     Znew[i][0] = sqrt(rsqnew) * cos(thetanew); */
/*     Znew[i][1] = sqrt(rsqnew) * sin(thetanew); */
/*     usesigma = Sigma[Ki[i] - 1]; */
    usesigma = zdelta;
    for(j=0;j<k;j++){ 
      Znew[i][j] = Z[i][j] + rnorm(0,sqrt(usesigma));
    }
    /*     Rprintf("Z[%d][%d] = %1.4f -->> Znew[%d][%d] = %1.4f\n",i,0,Z[i][0],i,0,Znew[i][0]);  */
  

    /* MOVE THIS TO NEW PROCRUSTES   -> SCALE ? */
    rms = 0.0;
    for(loop1=0;loop1<g;loop1++)
      for(j=0;j<k;j++)
	rms += Znew[loop1][j] * Znew[loop1][j];
    rms = sqrt(rms/((float)g*k));
    /*     Rprintf("RMS = %1.4f\n",rms); */
    /* ENDSCALE */

    for(loop1=0;loop1<g;loop1++)
      for(j=0;j<k;j++)
	Znew[loop1][j] = Znew[loop1][j]/rms;

    /*     Rprintf("About to procrustes transform Z\n"); */
    init_dmatrix(pZ,g,k,0.0); 
    if(procr_transform2(Znew,Z_mle,g,k,pZ,A, tZ, tZo, Ahalf, AhalfInv,
			tptrans, eAvectors, eADvalues, teAvectors, avZ,
			avZo, eAvalues, kk_helper, gg_helper, gk_helper,
			kg_helper, kk2_helper) == FOUND){
      init_dmatrix(Znew,g,k,0.0);
      copy_dmatrix(pZ,Znew,g,k);
    }
    /*     Rprintf("Procrustes transformed Z\n"); */

    /*     if(i < 5) */
    /*       Rprintf("mu[%d] = %1.4f %1.4f     sigma[%d] = %1.4f\n",Ki[i],mu[Ki[i] -1], mu[Ki[i]-1 + ng], Ki[i], Sigma[Ki[i]-1]); */
    llknew = loglike2_y(heads,tails,n_edges,g,Znew,k,beta,p,dir,X,mu,Sigma,Ki,ng);
    llkold = *llk;

    lr = llknew-llkold;
    /*     Rprintf("llkold = %1.4f    llknew = %1.4f    lr = %1.4f\n",llkold,llknew,lr); */
/*     loop2 = i; */
/*     sdlnorm(&longdim,ng,Ki[loop2]-1,mu,Sigma,Znew[loop2], &temp); */
/*     lr += temp; */
/*     sdlnorm(&longdim,ng,Ki[loop2]-1,mu,Sigma,Z[loop2], &temp); */
/*     lr -= temp; */

/*     Rprintf("Likelihood done\n"); */
    /*     Rprintf("Z[%d][%d] = %1.4f  <==> Znew[%d][%d] = %1.4f\n",i,0,Z[i][0],i,0,Znew[i][0]); */
    /*     Rprintf("Z[%d][%d] = %1.4f  <==> Znew[%d][%d] = %1.4f\n",i,1,Z[i][1],i,1,Znew[i][1]); */

    temp = runif(0.0,1.0); 
    if( temp < exp(lr) ){
      /*       Rprintf("Changing  ->%1.4f<-->%1.4f<-\n",temp,exp(lr)); */
      *llk = llknew;
      /*       for(j=0;j<k;j++) */
      /* 	Z[i][j] = Znew[i][j]; */
      copy_dmatrix(Znew,Z,g,k);
      change = 1;
    }
    else{
      /*       Rprintf("Not Changing\n"); */
      /*       for(j=0;j<k;j++) */
      /*       	Znew[i][j] = Z[i][j]; */
      copy_dmatrix(Z,Znew,g,k);
      *llk = llkold;
    }
    /*     Rprintf("Loop done\n"); */
  }/* end for i */

  if( change > 0 ){
    return(NEW);
  }
  else{
    return(OLD);
  }
  free_dmatrix(pZ,g);

/*   FREE procrustes stuff! */
/*   free_dmatrix(A,k); */
/*   free_dmatrix(tZ,k); */
/*   free_dmatrix(tZo,k); */
/*   free_dmatrix(gk_helper,g); */
/*   free_dmatrix(kk_helper,k); */
/*   free_dmatrix(gg_helper,g); */
/*   free_dmatrix(kk2_helper,k); */
/*   free_dmatrix(kg_helper,k); */
/*   free_dmatrix(eAvectors,k); */
/*   free_dmatrix(eADvalues,k); */
/*   free_dmatrix(teAvectors,k); */
/*   free_dmatrix(Ahalf,k); */
/*   free_dmatrix(AhalfInv,k); */
/*   free_dmatrix(tptrans,k); */
/*   free_dvector(avZ); */
/*   free_dvector(avZo); */
/*   free_dvector(eAvalues); */
/*   end free procrustes stuff! */
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
int beta_up2(int *heads, int *tails, int n_edges, double **Z, int g, int k, 
	     double prior_mu, double prior_sd, double *llk,
	     double *beta, int p, int dir, double ***X, double bdelta,
	     double *mu,double *Sigma,int *Ki,int ng)
{
  int i=0;
  double llknew=0.0, llkold=0.0, lr=0.0;
  double *betanew = dvector(p+1);

  for(i=0;i<=p;i++){
    betanew[i] = beta[i] + rnorm(0.0,bdelta);
  }

  llknew = loglike2_y(heads,tails,n_edges,g,Z,k,betanew,p,dir,X,mu,Sigma,Ki,ng);
  llkold = *llk;
  /* calculate the acceptance ratio */
  lr = llknew-llkold;
  /* this is the "Q" part */
  for(i=0;i<=p;i++){
    lr += dnorm(betanew[i],prior_mu,prior_sd,1);
    lr -= dnorm(beta[i],   prior_mu,prior_sd,1);  
  }
  /*   Rprintf("lr = %1.4f,  exp(lr) = %1.4f\n",lr,exp(lr)); */
  
  /* accept with probability r */
  if( runif(0.0,1.0) < exp(lr) ){
    for(i=0;i<=p;i++){
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
 CHANGED FROM Susan's CODE by passing in all necessary memory
*/
int procr_transform2(double **Z, double **Zo, int g, int k, double **pZ,
		     double **A, double **tZ, double **tZo, double ** Ahalf, 
		     double **AhalfInv, double **tptrans, double **eAvectors, 
		     double **eADvalues, double **teAvectors, double *avZ, 
		     double *avZo, double *eAvalues, double **kk_helper, 
		     double **gg_helper, double **gk_helper, 
		     double **kg_helper, double **kk2_helper)
{
  
  int i=0,j=0;
  
  /*   double *tempZ, *tempZo, *temp ; */

  /*  Center both Z and Zo around the mean of Zo    */
  /*  First centers Z around the origin  */
  
  for(j=0;j<k;j++){
    avZ[j] = 0.0;
    avZo[j] = 0.0;
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
  t(Z,g,k,tZ);
  t(Zo,g,k,tZo);

  init_dmatrix(A,k,k,0.0);
  init_dmatrix(kk_helper,k,k,0.0);
  init_dmatrix(kg_helper,k,g,0.0);
  dmatrix_multiply(tZ,k,g,Zo,k,kk_helper);
  dmatrix_multiply(kk_helper,k,k,tZo,g,kg_helper);
  dmatrix_multiply(kg_helper,k,g,Z,k,A);

  /* Compute sqrt(A) */
  sym_eigen(A,k,1,eAvalues,eAvectors);
  for(i=0;i<k;i++){
    eADvalues[i][i] = sqrt(eAvalues[i]);
  }
  t(eAvectors,k,k,teAvectors);  

  init_dmatrix(kk_helper,k,k,0.0);
  init_dmatrix(Ahalf,k,k,0.0);
  dmatrix_multiply(eAvectors,k,k,eADvalues,k,kk_helper);
  dmatrix_multiply(kk_helper,k,k,teAvectors,k,Ahalf);

  /*   init_dmatrix(AhalfInv,k,k,0.0); */
  /* Now compute the inverse */
  if(inverse(Ahalf,&k,AhalfInv) != FOUND) return NOTFOUND;

  /*Now compute t(Zo)*Z*AhalfInv*tZ*/
  init_dmatrix(kk_helper,k,k,0.0);
  init_dmatrix(kk2_helper,k,k,0.0);
  init_dmatrix(tptrans,k,g,0.0);
  dmatrix_multiply(tZo,k,g,Z,k,kk_helper);
  dmatrix_multiply(kk_helper,k,k,AhalfInv,k,kk2_helper);
  dmatrix_multiply(kk2_helper,k,k,tZ,g,tptrans);
  t(tptrans,k,g,pZ);

  return(FOUND);
}
