#ifndef CHANGESTATS_H
#define CHANGESTATS_H


#include "edgeTree.h"

/*  Node and dyad covariates are now passed as part of inputparams.  However, */
/*  attrib can still be set to point to the start of these attributes if */
/*  you want; see comments in InitErgmm.latent.r          Dave H  12/17/2003 */
struct OptionInput {
	void (*func)(int, Vertex*, Vertex*, struct OptionInput*, Gptr);
	double *attrib; /* Ptr to vector of attributes (node or dyad)*/
	int nstats;   /* Number of change statistics to be returned */
	double *dstats; /* ptr to change statistics returned */
	int ninputparams; /* Number of input parameters passed to function */
	double *inputparams; /* ptr to input parameters passed */
};


/*  change_statistics is an array of doubles of length n_param.
It is a temporary workspace for holding the value of the
change in the graph statistics from one graph to the next.
*/

#ifdef __BASECHANGESTATS_C__
#define extern
#endif
extern double *change_statistics;


void d_latentcov (int ntoggles, Vertex *heads, Vertex *tails, 
	      struct OptionInput *inp, Gptr g);

#endif
