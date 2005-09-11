#ifndef MODEL_H
#define MODEL_H

/*  Next to flag for Windows specific code */
/*  In reality this is done in ~/RCrossbuild/Makefile */
/*  under the ergmm entry. It appends these lines to the end. */
/*   */
/* #ifndef Win32 */
/*  #define Win32 */
/* #endif */

#include "edgeTree.h"
#include "globalStats.h"
#include "basechangeStats.h"

/* #ifndef OPTIONINPUT */
/* #define OPTIONINPUT */

/* struct OptionInput { */
/*   void (*func)(int, Vertex*, Vertex*, struct OptionInput*); */
/*   double *attrib;   Pointer to vector of attributes (node or dyad)*/
/*   int nstats;    Number of change statistics to be returned */
/*   double *dstats;  ptr to change statistics returned */
/*   int ninputparams;  Number of input parameters passed to function */
/*   double *inputparams;  ptr to input parameters passed */
/* };  */

/* #endif */

/* struct ModelOption is a lookup table tying a particular 
   model option (passed from R as a string) to a pointer to
   a function for computing change statistics in C.
*/
struct ModelOption {
  char *name;                                               /* option name */
  void (*func)(int, Vertex*, Vertex*, struct OptionInput*, Gptr g); /* function */
};

/* selected_options is modified to contain all necessary information
   for passing to various change-statistic-computing functions.
   Each such function (pointed to by selected_options.func[i])
   receives the entire OptionInput structure as part of its input
   and modifies the .dstats vector of statistics.
*/
extern struct OptionInput *selected_options;




/*  n_options is the number of OptionInput structures required to
    completely specify the model; total_options is the total number
    of different options available in option_list.
*/
extern int n_options;  
/* int total_options; Commented 12/15/03 by CTB */

/*  n_param is the total number of parameters in the model (which is
    the same as the number of graph statistics in the model).  Note that
    this is not necessarily the same as n_options, unless each option
    returns just a single statistic.
*/
extern int n_param;

int ModelInitialize (char *funnames, char *sonames, double *inputs);

void ModelDestroy();

void DegreeBoundInitialize(int *attribs, int *maxout, int *maxin, int *minout, int *minin, int condAllDegExact, int attriblength);

int GetIndexForAttrValue(int value);


#endif
