#ifndef PLINFO_H
#define PLINFO_H

#include "edgeTree.h"
#include "basechangeStats.h"
#include "model.h"

void plinfo_wrapper (double *heads, double *tails, double *dnedges,
                   double *dn, int *dflag, int *optionnum, char **funnames,
                   char **sonames, double *inputs,
                   double *responsevec, double *covmat,
                   int *start, int *end);

void plinfoInitialize (double *responsevec, double *covmat,
                       Vertex *start, Vertex *end,
                       Gptr g);

#endif
