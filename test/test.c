#include "hungarian.h"

#include <sys/time.h>
#include <time.h>

int main() {


  struct timeval tv1;
  struct timezone tzp1;  
  gettimeofday(&tv1, &tzp1);
  uint32_t random_seed = ((tv1.tv_sec & 0177) * 1000000) + tv1.tv_usec;

  fprintf(stderr, "random seed = %d\n", random_seed);
  srand(random_seed);
  
  uint32_t n = 100000;
  uint32_t m = 15;
  uint32_t i, j;
  
  WeightedBipartiteEdge *edges = malloc(n * m * sizeof(WeightedBipartiteEdge));

  uint32_t e = 0;
  for (i = 0; i < n; i++) {
    for(j = 0; j < m; j++) {
      edges[e].left = i;
      edges[e].right = rand() % n;
      edges[e].cost = j+1;
      e++;
    }
  }
  
  uint32_t *matching = hungarianMinimumWeightPerfectMatching(n, edges, n*m);

  if (matching == NULL) {
    fprintf(stdout, "Failure: Hungarian algorithm didn't find a matching.\n");
  } else {
    fprintf(stdout, "Matching was found.\n");
    for(i = 0; i < n; i++) {
      fprintf(stdout, "%u -> %u\n", i, matching[i]);
    }

    uint32_t weight = 0;
    for(i = 0; i < n*m; i++) {
      if(matching[edges[i].left] == edges[i].right) {
        weight += edges[i].cost;
      }
    }
    fprintf(stdout, "Minimum Weight is: %u\n", weight);
  }

  free(edges);
  free(matching);
  
  return 0;
}

