#ifndef GLOBALSTATS_H
#define GLOBALSTATS_H

#include "edgeTree.h"
#include "basechangeStats.h"

/* there are no current declarations for these functions....  */
/* but they being depricated anyway */

/* int undir_kstar (int k);  //will be replaced by sna calls?? */
/* int undir_triangles();//will be replaced by sna calls?? */

/*  do not allow any toggles that change the count of  */
/*  nodes above this degree. */
int g_condAboveDeg;

/* data to condition on degrees */
/*  condeg[degree][...]   lists of nodes with that degree */
Vertex **g_condeg;
/*  cached data to update condeg if we toggle */
int g_cOldHeadDeg; 
int g_cOldTailDeg; 
int g_cNewHeadDeg; 
int g_cNewTailDeg;
int g_fChange;

int g_fDistanceMetric;
int g_fMahalanobis;

/*  list of all pairs of nodes that have an edge between them. */
Vertex **g_rgAdjacent; /*  what edge did we replace?  */
int g_iEdge;
/*  total edges in graph */
int g_cEdgeCount;

/*  degcount[degree] = number of nodes with that degree */
int *g_degcount; 

int g_fInitCD;

int dist[65536];
int graphstate;

/*
  Variables for degree bounding
*/

int g_fBoundDegByAttr;
int g_attrcount;
int *g_attribs;
int *g_maxout;
int *g_maxin;
int *g_minout;
int *g_minin;

/*  condition on degree */
int OneRandomSwapCD (double *ratio);
int UpdateCDArrays();
int InitCDArrays();



double mahalanobis(double* graphstatistics, 
		   double* change_statistics, 
		   double* cholesky, 
		   int dimension);

double *cholesky;

void MCMC_global (double *heads, double *tails, double *dnedges,
		  double *dn, int *dflag, int *optionnum, char **funnames,
		  char **sonames, double *inputs,  
		  double *sample, 
		  int *fVerbose
		   );
#endif
