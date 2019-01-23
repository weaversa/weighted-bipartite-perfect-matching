#include "hungarian.h"

int main() {
  
  uint32_t n = 100;
  uint32_t m = 1000;
  uint32_t i;
  
  WeightedBipartiteEdge *edges = malloc(m * sizeof(WeightedBipartiteEdge));

  for (i = 0; i < m; i++) {
    edges[i].left = rand() % n;
    edges[i].right = rand() % n;
    edges[i].cost = (rand() % 100) + 1;
  }
  
  uint32_t *matching = hungarianMinimumWeightPerfectMatching(n, edges, m);

  free(edges);
  
  if (matching == NULL) {
    fprintf(stdout, "Failure: Hungarian algorithm didn't find a matching.\n");
  } else {
    fprintf(stdout, "Matching was found.\n");
  }

  free(matching);
  
  return 0;
}

