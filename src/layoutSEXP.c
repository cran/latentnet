/*
######################################################################
#
# layoutSEXP.c
#
# Written by Carter T. Butts <buttsc@uci.edu>
# Last Modified 11/11/06
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/network package
#
# This file contains routines related to layoutSEXP methods for network
# objects.
#
######################################################################
*/
 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include "utils.h"
#include "access.h"
#include "layoutSEXP.h"

/*TWO-DIMENSIONAL LAYOUT ROUTINES--------------------------------------*/

SEXP network_layout_fg_R(SEXP nw, SEXP niter, SEXP param, SEXP loc)
/*
Calculate a two-dimensional Fruchterman-Reingold layout for (symmetrized) 
adjacency matrix d.  Positions (stored in loc) should be initialized
prior to calling this routine.
*/
{
  double frk,t,ded,xd,yd;
  double *dx,*dy;
  double rf,af;
  int niteration;
  double maxdelta, volume, coolexp, repulserad; 
  long int n;
  int j,k;
  int i;
  
  /*Define various things*/
  n=networkSize(nw);
//  Rprintf("n %d\n",n);
  niteration=(int)REAL(niter)[0];
//  Rprintf("niteration %d\n",niteration);
  maxdelta=REAL(param)[0];
  volume=REAL(param)[1];
  coolexp=REAL(param)[2];
  repulserad=REAL(param)[3];
//  Rprintf("maxdelta %f %f %f %f\n",maxdelta,volume,coolexp,repulserad);

//  Rprintf("loc %f %f %f %f\n",REAL(loc)[0],REAL(loc)[1],REAL(loc)[n],REAL(loc)[n+1]);

  frk=sqrt(volume/(double)n); /*Define the F-R constant*/

  /*Allocate memory for transient structures*/
  dx=(double *)R_alloc(n,sizeof(double));
  dy=(double *)R_alloc(n,sizeof(double));
  /*Run the annealing loop*/
  for(i=niteration;i>=0;i--){
    /*Set the temperature (maximum move/iteration)*/
    t=maxdelta*pow(i/(double)niteration,coolexp);
    /*Clear the deltas*/
    for(j=0;j<n;j++){
      dx[j]=0.0;
      dy[j]=0.0;
    }
    /*Increment deltas for each undirected pair*/
    for(j=0;j<n;j++)
      for(k=j+1;k<n;k++){
        /*Obtain difference vector*/
        xd=REAL(loc)[  j]-REAL(loc)[  k];
        yd=REAL(loc)[n+j]-REAL(loc)[n+k];
        ded=sqrt(xd*xd+yd*yd);  /*Get dyadic euclidean distance*/
        xd/=ded;                /*Rescale differences to length 1*/
        yd/=ded;
        /*Calculate repulsive "force"*/
        rf=frk*frk*(1.0/ded-ded*ded/repulserad);
        dx[j]+=xd*rf;        /*Add to the position change vector*/
        dx[k]-=xd*rf;
        dy[j]+=yd*rf;
        dy[k]-=yd*rf;
        /*Calculate the attractive "force"*/
//  Rprintf("j %d k %d adj %d\n",j,k,isAdjacent(nw,j+1,k+1,0));
        if(isAdjacent(nw,j+1,k+1,0)||isAdjacent(nw,k+1,j+1,0)){
          af=ded*ded/frk;
          dx[j]-=xd*af;        /*Add to the position change vector*/
          dx[k]+=xd*af;
          dy[j]-=yd*af;
          dy[k]+=yd*af;
        }
      }
    /*Dampen motion, if needed, and move the points*/
    for(j=0;j<n;j++){
      ded=sqrt(dx[j]*dx[j]+dy[j]*dy[j]);
      if(ded>t){                 /*Dampen to t*/
        ded=t/ded;
        dx[j]*=ded;
        dy[j]*=ded;
      }
      REAL(loc)[  j]+=dx[j];               /*Update positions*/
      REAL(loc)[n+j]+=dy[j];
    }
//  Rprintf("%d of %d to go\n",i,niteration);
  }
//  Rprintf("loc %f %f %f %f\n",REAL(loc)[0],REAL(loc)[1],REAL(loc)[n],REAL(loc)[n+1]);
  return(loc);
}

int networkSize(SEXP x)
	/*Return the order of x.*/
{
	  SEXP atptr;

	    atptr = coerceVector(getNetworkAttribute(x,"n"),INTSXP);

	      return INTEGER(atptr)[0];
}

int isAdjacent(SEXP x, int vi, int vj, int naOmit)
/*Returns 0 if not adjacent, 1 if adjacent, and NA_INTEGER if missing and naOmit=0.*/
{
  SEXP mel,el,edge,endpts;
  int i,j,flag,isna,matchna,pc=0;
  
  /*Rprintf("\tInternal isAdjacent: seeking (%d,%d) w/naOmit=%d\n",vi,vj,naOmit);*/
  mel=getListElement(x,"mel");

  /*Start by hunting the outgoing edges of vi*/
  PROTECT(el=coerceVector(VECTOR_ELT(getListElement(x,"oel"),vi-1),INTSXP)); pc++;
  /*Rprintf("\t\tGot outgoing edge list for %d\n",vi);*/
  matchna=0;
  for(i=0;i<length(el);i++){
    /*Rprintf("\t\t\tCurrent edge is %d\n",INTEGER(el)[i]);*/
    edge=VECTOR_ELT(mel,INTEGER(el)[i]-1);
    isna=INTEGER(getListElement(getListElement(edge,"atl"),"na"))[0];  /*Missing edge?*/
    /*Rprintf("\t\t\tEdge missing status=%d\n",isna);*/
    PROTECT(endpts=coerceVector(getListElement(edge,"inl"),INTSXP)); pc++;      /*Get endpoints*/
    /*Rprintf("\t\t\tGot endpoints...looking for a match\n");*/
    flag=0;
    for(j=0;(!flag)&&(j<length(endpts));j++)   /*Check head of edge for vj*/
      if(INTEGER(endpts)[j]==vj){
        if(!isna){                             /*Return 1 on a clean match*/
          UNPROTECT(pc);
          return 1;
        }else{                  /*If matches but missing, note and move on*/
          matchna++;
          flag++;
        }
      }
  }
  /*Rprintf("\t\tDidn't find match...");*/
  /*Take stock of the situation*/
  if(isDirected(x)){            /*If directed, we're done here...*/
    if(matchna&&(!naOmit)){
      UNPROTECT(pc);
      return NA_INTEGER;          /*Matched a missing edge, not discounting*/
    }else{
      UNPROTECT(pc);
      return 0;                   /*Matched a missing edge, discounting*/
    }
  }

  /*If we're still here, x is undirected.  Now, we try vi's inedges.*/
  el=VECTOR_ELT(getListElement(x,"iel"),vi-1);
  PROTECT(el=coerceVector(VECTOR_ELT(getListElement(x,"iel"),vi-1),INTSXP)); pc++;
  for(i=0;i<length(el);i++){
    edge=VECTOR_ELT(mel,INTEGER(el)[i]-1);
    isna=INTEGER(getListElement(getListElement(edge,"atl"),"na"))[0];  /*Missing edge?*/
    PROTECT(endpts=coerceVector(getListElement(edge,"outl"),INTSXP)); pc++;      /*Get endpoints*/
    flag=0;
    for(j=0;(!flag)&&(j<length(endpts));j++)   /*Check tail of edge for vj*/
      if(INTEGER(endpts)[j]==vj){
        if(!isna){                             /*Return 1 on a clean match*/
          UNPROTECT(pc);
          return 1;
        }else{                  /*If matches but missing, note and move on*/
          matchna++;
          flag++;
        }
      }
  }
  /*Make the final decision*/
  if(matchna&&(!naOmit)){
    UNPROTECT(pc);
    return NA_INTEGER;          /*Matched a missing edge, not discounting*/
  }else{
    UNPROTECT(pc);
    return 0;                   /*Matched a missing edge, discounting*/
  }
}

int isDirected(SEXP x)
/*Is x directed?*/
{
  SEXP atptr;
  
  atptr = coerceVector(getNetworkAttribute(x,"directed"),LGLSXP);
  
  return INTEGER(atptr)[0];
}

SEXP getListElement(SEXP list, char *str)
/*Given a list and a character string, return a pointer to the element with the specified name (or else NULL).  This is taken from the Writing R Extensions manual.*/
{
  SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
  int i;
     
  for (i = 0; i < length(list); i++)
    if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
      elmt = VECTOR_ELT(list, i);
      break;
    }
  return elmt;
}

SEXP getNetworkAttribute(SEXP x, char *str)
/*Returns a pointer to the network attribute of x named by str, or else R_NilValue.*/
{
  return getListElement(getListElement(x,"gal"),str);
}
