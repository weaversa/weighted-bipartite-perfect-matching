// James Payor - December 2017
// MIT Licensed

// Rejiggered into C by Sean Weaver - January 2019

#include "hungarian.h"

// Macros!
//
// fo(i, n):       foreach i in [0, n)
// range(i, a, b): foreach i in [a, b)
//
// (Note: end value is cached, so fo(i, function()) will only have function called once.)

static const uint32_t oo = ~0;
static const uint32_t UNMATCHED = ~0;

typedef struct LeftEdge {
  uint32_t right;
  uint32_t cost;
} LeftEdge;

int LeftEdgeCmp (const void *a, const void *b) {
  LeftEdge arg1 = *(const LeftEdge*)a;
  LeftEdge arg2 = *(const LeftEdge*)b;

  if(arg1.right < arg2.right) return -1;
  if(arg1.right > arg2.right) return 1;
  if(arg1.cost < arg2.cost) return -1;
  if(arg1.cost > arg2.cost) return 1;
  return 0;
}

uint32_t *hungarianMinimumWeightPerfectMatching(uint32_t n, WeightedBipartiteEdge *allEdges, uint32_t numEdges) {

  // Edge lists for each left node.
  LeftEdge **leftEdges = malloc(n * sizeof(LeftEdge *));
  
  //region Edge list initialization
  
  // Initialize edge lists for each left node, based on the incoming set of edges.
  // While we're at it, we check that every node has at least one associated edge.
  // (Note: We filter out the edges that invalidly refer to a node on the left or right outside [0, n).)
  uint32_t *leftEdgeCounts = calloc(n, sizeof(uint32_t));
  uint32_t *rightEdgeCounts = calloc(n, sizeof(uint32_t));

  uint32_t edgeIndex;
  for(edgeIndex = 0; edgeIndex < numEdges; edgeIndex++) {
    WeightedBipartiteEdge edge = allEdges[edgeIndex];
    if (edge.left >= 0 && edge.left < n) {
      leftEdgeCounts[edge.left]++;
    }
    if (edge.right >= 0 && edge.right < n) {
      rightEdgeCounts[edge.right]++;
    }
  }

  uint32_t i;
  for(i = 0; i < n; i++) {
    if (leftEdgeCounts[i] == 0 || rightEdgeCounts[i] == 0) {
      // No matching will be possible.

      free(leftEdgeCounts);
      free(rightEdgeCounts);
      free(leftEdges);
      
      return NULL;
    }
  }

  free(rightEdgeCounts);
  
  // Reserve the required space for each node.
  for(i = 0; i < n; i++) {
    leftEdges[i] = malloc(leftEdgeCounts[i] * sizeof(LeftEdge));
    leftEdgeCounts[i] = 0;
  }
  // Actually add to the edge lists now.
  for(edgeIndex = 0; edgeIndex < numEdges; edgeIndex++) {
    WeightedBipartiteEdge edge = allEdges[edgeIndex];
    if (edge.left >= 0 && edge.left < n && edge.right >= 0 && edge.right < n) {
      leftEdges[edge.left][leftEdgeCounts[edge.left]++] = (struct LeftEdge) {edge.right, edge.cost};
    }
  }

  // Sort the edge lists, and remove duplicate edges (keep the edge with smallest cost).
  for(i = 0; i < n; i++) {
    LeftEdge *edges = leftEdges[i];

    qsort(edges, leftEdgeCounts[i], sizeof(LeftEdge), LeftEdgeCmp);
    
    int edgeCount = 0;
    for(edgeIndex = 1; edgeIndex < leftEdgeCounts[i]; edgeIndex++) {
      LeftEdge edge = edges[edgeIndex];
      if (edge.right != edges[edgeCount].right) {
	edges[++edgeCount] = edge;
      }
    }

    leftEdgeCounts[i] = edgeCount+1;
  }

  //endregion Edge list initialization
  
  // These hold "potentials" for nodes on the left and nodes on the
  // right, which reduce the costs of attached edges.  We maintain
  // that every reduced cost, cost[i][j] - leftPotential[i] -
  // leftPotential[j], is greater than zero.
  uint32_t *leftPotential = malloc(n * sizeof(uint32_t));
  uint32_t *rightPotential = malloc(n * sizeof(uint32_t));
  
  //region Node potential initialization
  
  // Here, we seek good initial values for the node potentials.
  // Note: We're guaranteed by the above code that at every node on the left and right has at least one edge.
  
  // First, we raise the potentials on the left as high as we can for each node.
  // This guarantees each node on the left has at least one "tight" edge.
  
  for(i = 0; i < n; i++) {
    LeftEdge *edges = leftEdges[i];
    uint32_t smallestEdgeCost = edges[0].cost;
    for(edgeIndex = 1; edgeIndex < leftEdgeCounts[i]; edgeIndex++) {
      if (edges[edgeIndex].cost < smallestEdgeCost) {
	smallestEdgeCost = edges[edgeIndex].cost;
      }
    }
    
    // Set node potential to the smallest incident edge cost.
    // This is as high as we can take it without creating an edge with zero reduced cost.
    leftPotential[i] = smallestEdgeCost;
  }
  
  // Second, we raise the potentials on the right as high as we can for each node.
  // We do the same as with the left, but this time take into account that costs are reduced
  // by the left potentials.
  // This guarantees that each node on the right has at least one "tight" edge.

  memset(rightPotential, 0xff, n * sizeof(uint32_t));

  for(i = 0; i < n; i++) {
    LeftEdge *edges = leftEdges[i];
    for(edgeIndex = 0; edgeIndex < leftEdgeCounts[i]; edgeIndex++) {
      LeftEdge edge = edges[edgeIndex];
      uint32_t reducedCost = edge.cost - leftPotential[i];
      assert(reducedCost <= edge.cost);
      if (rightPotential[edge.right] > reducedCost) {
	rightPotential[edge.right] = reducedCost;
      }
    }
  }

  //endregion Node potential initialization

  // Tracks how many edges for each left node are "tight".
  // Following initialization, we maintain the invariant that these are at the start of the node's edge list.
  uint32_t *leftTightEdgesCount = calloc(n, sizeof(uint32_t));
  
  //region Tight edge initialization
  
  // Here we find all tight edges, defined as edges that have zero reduced cost.
  // We will be interested in the subgraph induced by the tight edges, so we partition the edge lists for
  // each left node accordingly, moving the tight edges to the start.

  for(i = 0; i < n; i++) {
    LeftEdge *edges = leftEdges[i];
    uint32_t tightEdgeCount = 0;
    for(edgeIndex = 0; edgeIndex < leftEdgeCounts[i]; edgeIndex++) {
      LeftEdge edge = edges[edgeIndex];
      uint32_t reducedCost = edge.cost - leftPotential[i] - rightPotential[edge.right];
      assert(reducedCost <= edge.cost);
      if (reducedCost == 0) {
	if (edgeIndex != tightEdgeCount) {
	  //Swap edges
	  LeftEdge tmp = edges[tightEdgeCount];
	  edges[tightEdgeCount] = edges[edgeIndex];
	  edges[edgeIndex] = tmp;
	}
	tightEdgeCount++;
      }
    }
    leftTightEdgesCount[i] = tightEdgeCount;
  }
  
  //endregion Tight edge initialization

  // Now we're ready to begin the inner loop.
  
  // We maintain an (initially empty) partial matching, in the subgraph of tight edges.
  uint32_t currentMatchingCardinality = 0;

  uint32_t *leftMatchedTo = malloc(n * sizeof(uint32_t));
  uint32_t *rightMatchedTo = malloc(n * sizeof(uint32_t));
  memset(leftMatchedTo, 0xff, n * sizeof(uint32_t));
  memset(rightMatchedTo, 0xff, n * sizeof(uint32_t));
  
  //region Initial matching (speedup?)
  
  // Because we can, let's make all the trivial matches we can.
  for(i = 0; i < n; i++) {
    LeftEdge *edges = leftEdges[i];
    for(edgeIndex = 0; edgeIndex < leftTightEdgesCount[i]; edgeIndex++) {
      uint32_t j = edges[edgeIndex].right;
      if (rightMatchedTo[j] == UNMATCHED) {
	currentMatchingCardinality++;
	rightMatchedTo[j] = i;
	leftMatchedTo[i] = j;
	break;
      }
    }
  }
  
  if (currentMatchingCardinality == n) {
    // Well, that's embarassing. We're already done!
    free(leftEdgeCounts);

    for(i = 0; i < n; i++) {
      free(leftEdges[i]);
    }
    free(leftEdges);

    free(leftPotential);
    free(rightPotential);
    
    free(rightMatchedTo);
      
    return leftMatchedTo;
  }

  return NULL;
  
  /*
  
  //endregion Initial matching (speedup?)
  
  // While an augmenting path exists, we add it to the matching.
  // When an augmenting path doesn't exist, we update the potentials so that an edge between the area
  // we can reach and the unreachable nodes on the right becomes tight, giving us another edge to explore.
  //
  // We proceed in this fashion until we can't find more augmenting paths or add edges.
  // At that point, we either have a min-weight perfect matching, or no matching exists.
  
  //region Inner loop state variables
  
  // One point of confusion is that we're going to cache the edges between the area we've explored
  // that are "almost tight", or rather are the closest to being tight.
  // This is necessary to achieve our O(N^3) runtime.
  //
  // rightMinimumSlack[j] gives the smallest amount of "slack" for an unreached node j on the right,
  // considering the edges between j and some node on the left in our explored area.
  //
  // rightMinimumSlackLeftNode[j] gives the node i with the corresponding edge.
  // rightMinimumSlackEdgeIndex[j] gives the edge index for node i.
  
  int rightMinimumSlack[n], rightMinimumSlackLeftNode[n], rightMinimumSlackEdgeIndex[n];
  
  std::deque<int> leftNodeQueue;
  bool leftSeen[n];
  int rightBacktrack[n];
  
  // Note: the above are all initialized at the start of the loop.
  
  //endregion Inner loop state variables
  
  while (currentMatchingCardinality < n) {
    
    //region Loop state initialization
    
    // Clear out slack caches.
    // Note: We need to clear the nodes so that we can notice when there aren't any edges available.
    std::fill_n(rightMinimumSlack, n, oo);
    std::fill_n(rightMinimumSlackLeftNode, n, UNMATCHED);
    
    // Clear the queue.
    leftNodeQueue.clear();
    
    // Mark everything "unseen".
    std::fill_n(leftSeen, n, false);
    std::fill_n(rightBacktrack, n, UNMATCHED);
    
    //endregion Loop state initialization
    
    int startingLeftNode = UNMATCHED;
    
    //region Find unmatched starting node
    
    // Find an unmatched left node to search outward from.
    // By heuristic, we pick the node with fewest tight edges, giving the BFS an easier time.
    // (The asymptotics don't care about this, but maybe it helps. Eh.)
    {
      int minimumTightEdges = oo;
      fo(i, n) {
	if (leftMatchedTo[i] == UNMATCHED && leftTightEdgesCount[i] < minimumTightEdges) {
	  minimumTightEdges = leftTightEdgesCount[i];
	  startingLeftNode = i;
	}
      }
    }
    
    //endregion Find unmatched starting node
    
    assert(startingLeftNode != UNMATCHED);
    
    assert(leftNodeQueue.empty());
    leftNodeQueue.push_back(startingLeftNode);
    leftSeen[startingLeftNode] = true;
    
    int endingRightNode = UNMATCHED;
    while (endingRightNode == UNMATCHED) {
      
      //region BFS until match found or no edges to follow
      
      while (endingRightNode == UNMATCHED && !leftNodeQueue.empty()) {
	// Implementation note: this could just as easily be a DFS, but a BFS probably
	// has less edge flipping (by my guess), so we're using a BFS.
	
	const int i = leftNodeQueue.front();
	leftNodeQueue.pop_front();
	
	std::vector<LeftEdge>& edges = leftEdges[i];
	// Note: Some of the edges might not be tight anymore, hence the awful loop.
	for(int edgeIndex = 0; edgeIndex < leftTightEdgesCount[i]; ++edgeIndex) {
	  const LeftEdge& edge = edges[edgeIndex];
	  const int j = edge.right;
	  
	  assert(edge.cost - leftPotential[i] - rightPotential[j] >= 0);
	  if (edge.cost > leftPotential[i] + rightPotential[j]) {
	    // This edge is loose now.
	    --leftTightEdgesCount[i];
	    std::swap(edges[edgeIndex], edges[leftTightEdgesCount[i]]);
	    --edgeIndex;
	    continue;
	  }
	  
	  if (rightBacktrack[j] != UNMATCHED) {
	    continue;
	  }
	  
	  rightBacktrack[j] = i;
	  int matchedTo = rightMatchedTo[j];
	  if (matchedTo == UNMATCHED) {
	    // Match found. This will terminate the loop.
	    endingRightNode = j;
	    
	  } else if (!leftSeen[matchedTo]) {
	    // No match found, but a new left node is reachable. Track how we got here and extend BFS queue.
	    leftSeen[matchedTo] = true;
	    leftNodeQueue.push_back(matchedTo);
	  }
	}
	
	//region Update cached slack values
	
	// The remaining edges may be to nodes that are unreachable.
	// We accordingly update the minimum slackness for nodes on the right.
	
	if (endingRightNode == UNMATCHED) {
	  const int potential = leftPotential[i];
	  range(edgeIndex, leftTightEdgesCount[i], edges.size()) {
	    const LeftEdge& edge = edges[edgeIndex];
	    int j = edge.right;
	    
	    if (rightMatchedTo[j] == UNMATCHED || !leftSeen[rightMatchedTo[j]]) {
	      // This edge is to a node on the right that we haven't reached yet.
	      
	      int reducedCost = edge.cost - potential - rightPotential[j];
	      assert(reducedCost >= 0);
	      
	      if (reducedCost < rightMinimumSlack[j]) {
		// There should be a better way to do this backtracking...
		// One array instead of 3. But I can't think of something else. So it goes.
		rightMinimumSlack[j] = reducedCost;
		rightMinimumSlackLeftNode[j] = i;
		rightMinimumSlackEdgeIndex[j] = edgeIndex;
	      }
	    }
	  }
	}
	
	//endregion Update cached slack values
      }
      
      //endregion BFS until match found or no edges to follow
      
      //region Update node potentials to add edges, if no match found
      
      if (endingRightNode == UNMATCHED) {
	// Out of nodes. Time to update some potentials.
	int minimumSlackRightNode = UNMATCHED;
	
	//region Find minimum slack node, or abort if none exists
	
	int minimumSlack = oo;
	fo(j, n) {
	  if (rightMatchedTo[j] == UNMATCHED || !leftSeen[rightMatchedTo[j]]) {
	    // This isn't a node reached by our BFS. Update minimum slack.
	    if (rightMinimumSlack[j] < minimumSlack) {
	      minimumSlack = rightMinimumSlack[j];
	      minimumSlackRightNode = j;
	    }
	  }
	}
	
	if (minimumSlackRightNode == UNMATCHED || rightMinimumSlackLeftNode[minimumSlackRightNode] == UNMATCHED) {
	  // The caches are all empty. There was no option available.
	  // This means that the node the BFS started at, which is an unmatched left node, cannot reach the
	  // right - i.e. it will be impossible to find a perfect matching.
	  
	  return std::vector<int>();
	}
	
	//endregion Find minimum slack node, or abort if none exists
	
	assert(minimumSlackRightNode != UNMATCHED);
	
	// Adjust potentials on left and right.
	fo(i, n) {
	  if (leftSeen[i]) {
	    leftPotential[i] += minimumSlack;
	    if (leftMatchedTo[i] != UNMATCHED) {
	      rightPotential[leftMatchedTo[i]] -= minimumSlack;
	    }
	  }
	}
	
	// Downward-adjust slackness caches.
	fo(j, n) {
	  if (rightMatchedTo[j] == UNMATCHED || !leftSeen[rightMatchedTo[j]]) {
	    rightMinimumSlack[j] -= minimumSlack;
	    
	    // If the slack hit zero, then we just found ourselves a new tight edge.
	    if (rightMinimumSlack[j] == 0) {
	      const int i = rightMinimumSlackLeftNode[j];
	      const int edgeIndex = rightMinimumSlackEdgeIndex[j];
	      
	      //region Update leftEdges[i] and leftTightEdgesCount[i]
	      
	      // Move it in the relevant edge list.
	      if (edgeIndex != leftTightEdgesCount[i]) {
		std::vector<LeftEdge>& edges = leftEdges[i];
		std::swap(edges[edgeIndex], edges[leftTightEdgesCount[i]]);
	      }
	      ++leftTightEdgesCount[i];
	      
	      //endregion Update leftEdges[i] and leftTightEdgesCount[i]
	      
	      // If we haven't already encountered a match, we follow the edge and update the BFS queue.
	      // It's possible this edge leads to a match. If so, we'll carry on updating the tight edges,
	      // but won't follow them.
	      if (endingRightNode == UNMATCHED) {
		// We're contemplating the consequences of following (i, j), as we do in the BFS above.
		rightBacktrack[j] = i;
		int matchedTo = rightMatchedTo[j];
		if (matchedTo == UNMATCHED) {
		  // Match found!
		  endingRightNode = j;
		} else if (!leftSeen[matchedTo]) {
		  // No match, but new left node found. Extend BFS queue.
		  leftSeen[matchedTo] = true;
		  leftNodeQueue.push_back(matchedTo);
		}
	      }
	    }
	  }
	}
      }
      
      //endregion Update node potentials to add edges, if no match found
    }
    
    // At this point, we've found an augmenting path between startingLeftNode and endingRightNode.
    // We'll just use the backtracking info to update our match information.
    
    ++currentMatchingCardinality;
    
    //region Backtrack and flip augmenting path
    
    {
      int currentRightNode = endingRightNode;
      while (currentRightNode != UNMATCHED) {
	const int currentLeftNode = rightBacktrack[currentRightNode];
	const int nextRightNode = leftMatchedTo[currentLeftNode];
	
	rightMatchedTo[currentRightNode] = currentLeftNode;
	leftMatchedTo[currentLeftNode] = currentRightNode;
	
	currentRightNode = nextRightNode;
      }
    }
    
    //endregion Backtrack and flip augmenting path
  }
  
  // Oh look, we're done.
  return std::vector<int>(leftMatchedTo, leftMatchedTo + n);
  */
}
