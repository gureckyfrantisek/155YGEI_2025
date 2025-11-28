from lines_to_graph2 import *
from pathlib import Path
from r2g_ver_1_2 import get_graphs, ADRESAR

def dijkstra(G, start):
    '''
        Dijkstra's algorithm
        G - graph on which we are operating
        start - the starting node
        
        returns a list of predecessors
    '''
    
    # Init variables
    n = len(G)
    pred = [None] * (n + 1)
    dist = [inf] * (n + 1)
    prioQ = PriorityQueue()
    
    # Put the start into queue and update the distance
    prioQ.put((0, start))
    dist[start] = 0
    
    # While we have nodes
    while not prioQ.empty():
        # Get the top prio node
        du, u = prioQ.get()
        
        # For all it's neighbours
        for v, wuv in G[u].items():
            # Relax the node if there's a shorter path
            if dist[v] > dist[u] + wuv:
                dist[v] = dist[u] + wuv
                pred[v] = u
                prioQ.put((dist[v], v))
                
    return pred

def rectPath(pred, u, v):
    path = []       # Empty path
    
    while u != v and v != None:
        path.append(v)
        v = pred[v]
        
    path.append(v)
    return path

def main():
    # Get graphs from OSM
    get_graphs()

    # Translate all edge data to graphs
    data_dir = Path(__file__).parent / ADRESAR
    graphs = []

    for file in data_dir.iterdir():
        if file.name == "nodes.txt":
            continue
        graphs.append(convertToGraphWithIDs(file))
    
    # Run Dijsktra on all costs
    predecessors = [[]] * len(graphs)
    paths = [[]] * len(graphs)
    
    for i in range(len(graphs)):
        predecessors[i] = dijkstra(graphs[i], 1)
        
        paths[i] = rectPath(predecessors[i], 1, 57851)
    
        print(paths[i])

if __name__ == "__main__":
    main()
