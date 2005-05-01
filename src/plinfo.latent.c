#include "plinfo.latent.h"

void plinfo_wrapper (double *heads, double *tails, double *dnedges,
		   double *dn, int *dflag, int *optionnum, char **funnames,
		   char **sonames, double *inputs,  
		   double *responsevec, double *covmat,
                   int *start, int *end)
{
  Gptr g;

  /*  need to ensure this is true - this flag is not always reset  */
  /*  (for instance, in case of errors that cause aborts midprogram) */
  g_fDistanceMetric = 0;

  GetRNGstate(); /* Necessary for setting up use of the R random number generator */
  n_nodes = (Vertex)*dn; /* Coerce double *dn to type Vertex */
  n_options = *optionnum;
  directed_flag = *dflag;

  g=GraphInitialize(heads, tails, (Edge)*dnedges);  /* Coerce *dnedges to type Edge */
  ModelInitialize(*funnames, *sonames, inputs);
  
  plinfoInitialize(responsevec, covmat, (Vertex*)start, (Vertex*)end, g);
  
  ModelDestroy();
  GraphDestroy();
  
  PutRNGstate(); /* Must be called after GetRNGstate before returning to R */
}

void plinfoInitialize (double *responsevec, double *covmat,
                       Vertex *start, Vertex *end,
                       Gptr g)
{
  int k, l, outflag, inflag = 0;
  Edge offset1 = (Edge) (n_nodes * (n_nodes-1))/2;
  Edge offset2 = (Edge) offset1 * n_param;
  Vertex i, j, dc;
  struct OptionInput *inp;

  dc = 0;
  for(i=1; i<n_nodes; i++){
   for(j = i+1; j <= n_nodes; j++){
    dc++;
    if((dc >= (*start)) & (dc <= (*end))){
     if (*(g.directed_flag)){
       *(responsevec + offset1) = inflag = (EdgetreeSearch(i, j, g.inedges) != 0);
     }
     *responsevec++ = outflag = (EdgetreeSearch(i, j, g.outedges) != 0);
     
     for (k=0, inp = selected_options; k < n_options; k++, inp++)
      {
        inp->dstats = covmat;
        (*(inp->func))(1, &i, &j, inp, g);
        
        /* dstats are returned for TOGGLED edges; for MPLE, we need values
           reflecting a change from 0->1 instead.  Thus, we have to change 
           the sign of dstats if the edge exists. */
        
        if (outflag)
         {
          for (l = 0; l < inp->nstats; l++)
           inp->dstats[l] = -inp->dstats[l];
         }
        
        if (*(g.directed_flag))
         {
         /*  then get the graph stats for these stats */
         inp->dstats = covmat + offset2;
         (*(inp->func))(1, &j, &i, inp, g);
         if (inflag)
          {
           for (l=0; l<inp->nstats; l++)
            inp->dstats[l] = -inp->dstats[l];
          }
         }
        covmat += inp->nstats;    
      }
    }
   }
  }
}
