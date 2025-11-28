from lines_to_graph2 import *
from pathlib import Path
from r2g_ver_1_2 import get_graphs, ADRESAR
import matplotlib.pyplot as plt
from kruskal import kruskal_mst

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

def show_graph(graph, paths=[]):
    # Prepare the dictionary
    coords = {}

    # Open and load the nodes.txt file with coordinates
    with open('data_grafy/nodes.txt') as f:
        for line in f:
            # Split with whitespace
            id, x, y = line.split()

            # Create the point in dictionary
            coords[int(id)] = [float(x), float(y)]
    
    # Show the graph as it is, plot every edge
    for node in graph:
        for neighbour in graph[node]:
            plt.plot([coords[node][0], coords[neighbour][0]], [coords[node][1], coords[neighbour][1]], color="blue", linewidth=1)

    # Highlight the paths for different costs
    for i, path in enumerate(paths):
        path_coords = [[], []]
        for node in path:
            path_coords[0].append(coords[node][0])
            path_coords[1].append(coords[node][1])

        plt.plot(path_coords[0], path_coords[1], linewidth=3, label=i)

    # Finally show the plot
    plt.axis('equal')
    plt.legend()
    plt.show()
            

def main():
    # Get graphs from OSM
    get_graphs()

    # Translate all edge data to graphs
    data_dir = Path(__file__).parent / ADRESAR
    graphs = []

    # Bonus: prepare arrays for minimum spanning trees
    kruskal_msts = []

    for file in data_dir.iterdir():
        if file.name == "nodes.txt":
            continue
        graphs.append(convertToGraphWithIDs(file))
        kruskal_msts.append(kruskal_mst(file))

    # Run Dijsktra on all costs
    predecessors = [[]] * len(graphs)
    paths = [[]] * len(graphs)
    
    for i in range(len(graphs)):
        predecessors[i] = dijkstra(graphs[i], 1)
        
        paths[i] = rectPath(predecessors[i], 1, 33456)
    
        print(paths[i])
    
    # Show the best paths for all weights
    # show_graph(graphs[0], paths)

    # Show the MSTs
    print(kruskal_msts[0])
    for mst in kruskal_msts:
        show_graph(mst)

if __name__ == "__main__":
    main()
