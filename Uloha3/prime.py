from queue import PriorityQueue
from collections import *

def prime_mst(G, start):
    # Init variables
    visited = set()
    prioQ = PriorityQueue()
    MST = defaultdict(dict)

    # Put the start into queue
    prioQ.put((0, start, None))  # (weight, current_node, parent_node)
    
    # While we have nodes
    while not prioQ.empty():
        # Get the top prio edge
        w, u, parent = prioQ.get()
        
        # Skip if already visited
        if u in visited:
            continue
        
        visited.add(u)
        
        # Add edge to MST if not the start node
        if parent is not None:
            MST[u][parent] = w
            MST[parent][u] = w
        
        # For all its neighbours
        for v, wuv in G[u].items():
            # Add unvisited neighbours to queue
            if v not in visited:
                prioQ.put((wuv, v, u))
    
    return MST