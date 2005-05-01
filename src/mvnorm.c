/****************************************************************************/
/*  Author: Jeremy Tantrum, tantrum@stat.washington.edu                     */
/*  Purpose: support functions for model 2                                  */
/*           proposed by Adrian E. Raftery, Mark S. Handcock and Jeremy T   */
/****************************************************************************/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <R.h>
#include <Rmath.h>

#include "latentUtil.h"


void sdlnorm(long *p, int ng, int grp, double *mu, double *Sigma, double *x, 
	     double *result)
{
  /*   -0.5[ (x-mu)^t Sig^-1 (x-mu) + *p * log(2 * Pi * det(Sig)) ] */
  int i;

  double jvarstore1;

  jvarstore1 = 0.0;
  for(i=0;i<(*p);i++)
    jvarstore1 += (x[i] - mu[grp + i * ng])*(x[i] - mu[grp + i * ng]);
  jvarstore1 = jvarstore1 / (-2 * Sigma[grp]);

  *result = -0.5 * (*p) * log(2 * M_PI * Sigma[grp]) + jvarstore1;
}


void mvdlnorm2(long *p, double *mu, double *Sigma, double *x, double *result);

double call_mvdlnorm2(long p, double *mu, double *Sigma, double *x)
{
  double result;
  mvdlnorm2(&p, mu, Sigma, x, &result);

  return(result);
}

void mvdlnorm2(long *p, double *mu, double *Sigma, double *x, double *result)
{
  /*   -0.5[ (x-mu)^t Sig^-1 (x-mu) + *p * log(2 * Pi * det(Sig)) ] */
  double vec1a, vec1b, vec2a, vec2b, si1,si2,si4, s1, s2, s4, 
    jvarstore1, jvarstore2;

  s1 = Sigma[0];
  s2 = Sigma[1];
  s4 = Sigma[3];

  jvarstore1 = s1 * s4 - s2 * s2;
  /*             Rprintf("det Sigma = %1.4f\n",jvarstore1); */
  si1 = s4/jvarstore1;
  /*             Rprintf("SigInv[0][0] = %1.4f\n",Siginv[0][0]); */
  si4 = s1/jvarstore1;
  /*             Rprintf("SigInv[0][0] = %1.4f\n",Siginv[1][1]); */
  si2 = -s2/jvarstore1;
  /*             Rprintf("SigInv[0][0] = %1.4f\n",Siginv[1][0]); */

  /*   Rprintf("det Sigma = %1.8g \n", jvarstore1); */

  jvarstore2 = -log(jvarstore1);

  vec1a = x[0] - mu[0];
  vec1b = x[1] - mu[1];

  /*   Rprintf("x-mu: (%1.4g, %1.4g)\n",vec1a,vec1b); */

  vec2a = vec1a * si1 + vec1b * si2;
  vec2b = vec1a * si2 + vec1b * si4;

  /*   Rprintf("(x-mu) * Siginv: (%1.4g, %1.4g)\n",vec2a,vec2b); */

  jvarstore1 = vec1a*vec2a + vec1b*vec2b;

  /*   Rprintf("log det = %1.4f || log of exp = %1.4f\n",jvarstore2, jvarstore1); */

  *result = -0.5 * ( jvarstore1 + (*p) * (jvarstore2 + log(2 * M_PI)));

  /*   Rprintf("result = %1.4f\n",*result); */
}
