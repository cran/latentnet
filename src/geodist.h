#ifndef GEODIST_H
#define GEODIST_H

#include "edgeTree.h"

void node_geodesics (Vertex *edgelist, Vertex *nnodes, Edge *nodelist,
                     Edge *nedges, int *nodecolor, Vertex *dist, 
                     Vertex *Q, Vertex *source);

void geodesic_matrix (Vertex *edgelist, Vertex *nnodes,
		      Edge *nodelist, Edge *nedges, 
		      int *nodecolor, Vertex *distmat, Vertex *Q);

#endif
