#ifndef GLOBALSTATS_H
#define GLOBALSTATS_H


#include "edgeTree.h"
#include "basechangeStats.h"

#ifdef __MODEL_C__
#define extern 
#endif

/* there are no current declarations for these functions....  */
/* but they being depricated anyway */

/* int undir_kstar (int k);  //will be replaced by sna calls?? */
/* int undir_triangles();//will be replaced by sna calls?? */

/*  do not allow any toggles that change the count of  */
/*  nodes above this degree. */
extern int g_condAboveDeg;

/* data to condition on degrees */
/*  condeg[degree][...]   lists of nodes with that degree */
extern Vertex **g_condeg;
/*  cached data to update condeg if we toggle */
extern int g_cOldHeadDeg; 
extern int g_cOldTailDeg; 
extern int g_cNewHeadDeg; 
extern int g_cNewTailDeg;
extern int g_fChange;

extern int g_fDistanceMetric;
extern int g_fMahalanobis;

/*  list of all pairs of nodes that have an edge between them. */
extern Vertex **g_rgAdjacent; /*  what edge did we replace?  */
extern int g_iEdge;
/*  total edges in graph */
extern int g_cEdgeCount;

/*  degcount[degree] = number of nodes with that degree */
extern int *g_degcount; 

extern int g_fInitCD;

extern int dist[65536];
extern int graphstate;

/*
  Variables for degree bounding
*/

extern int g_fBoundDegByAttr;
extern int g_attrcount;
extern int *g_attribs;
extern int *g_maxout;
extern int *g_maxin;
extern int *g_minout;
extern int *g_minin;

/*  condition on degree */
extern int OneRandomSwapCD (double *ratio);
extern int UpdateCDArrays();
extern int InitCDArrays();



double mahalanobis(double* graphstatistics, 
		   double* change_statistics, 
		   double* cholesky, 
		   int dimension);

extern double *cholesky;

void MCMC_global (double *heads, double *tails, double *dnedges,
		  double *dn, int *dflag, int *optionnum, char **funnames,
		  char **sonames, double *inputs,  
		  double *sample, 
		  int *fVerbose
		   );
#endif
