#ifndef HUNGARIAN_H
#define HUNGARIAN_H

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>

typedef struct WeightedBipartiteEdge {
  uint32_t left;
  uint32_t right;
  uint32_t cost;
} WeightedBipartiteEdge;

// Given the number of nodes on each side of the bipartite graph and a list of edges, returns a minimum-weight perfect matching.
// If a matching is found, returns a length-n vector, giving the nodes on the right that the left nodes are matched to.
// If no matching exists, returns an empty vector.
// (Note: Edges with endpoints out of the range [0, n) are ignored.)
uint32_t *hungarianMinimumWeightPerfectMatching(uint32_t n, WeightedBipartiteEdge *allEdges, uint32_t numEdges);

#endif
