#define __BASECHANGESTATS_C__

#include "basechangeStats.h"

void d_latentcov (int ntoggles, Vertex *heads, Vertex *tails, 
	      struct OptionInput *inp, Gptr g) 
{
  double val;
  Vertex h, t;
  int i, edgeflag, refedgeflag;

  for(i=0;i<3;i++)
    inp->dstats[i] = 0.0;
  for (i=0; i<ntoggles; i++) 
    {
      /*Get the initial state of the edge and its reflection*/
      edgeflag=(EdgetreeSearch(h=heads[i], t=tails[i], g.outedges) != 0);
      refedgeflag = (EdgetreeSearch(t, h, g.outedges) != 0);

      /*Get the dyadic covariate*/
      val = inp->inputparams[1+(t-1)+(h-1)*((long int)(inp->inputparams[0]))];
      
      /*Update the change statistics, as appropriate*/
      if(refedgeflag){      /* Reflected edge is present */
        if(edgeflag){         /* Toggled edge _was_ present */
	  if(t>h){              /* Mut to low->high */
	    inp->dstats[0] -= val;
	    inp->dstats[1] += val;
	  }else{                /* Mut to high->low */
	    inp->dstats[0] -= val;
	    inp->dstats[2] += val;
	  }
        }else{                /* Toggled edge _was not_ present */
	  if(t>h){              /* Low->high to mut */
	    inp->dstats[1] -= val;
	    inp->dstats[0] += val;
	  }else{                /* High->low to mut */
	    inp->dstats[2] -= val;
	    inp->dstats[0] += val;
	  }
	}
      }else{                /* Reflected edge is absent */
        if(edgeflag){         /* Toggled edge _was_ present */
	  if(t>h){              /* High->low to null */
	    inp->dstats[2] -= val;
	  }else{                /* Low->high to null */
	    inp->dstats[1] -= val;
	  }
        }else{                /* Toggled edge _was not_ present */
	  if(t>h){              /* Null to high->low */
	    inp->dstats[2] += val;
	  }else{                /* Null to low->high */
	    inp->dstats[1] += val;
	  }
	}
      }
      
      if (i+1 < ntoggles)
	ToggleEdge(heads[i], tails[i], g);  /* Toggle this edge if more to come */
    }
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    ToggleEdge(heads[i], tails[i], g); 
}

#undef __BASECHANGESTATS_C__
