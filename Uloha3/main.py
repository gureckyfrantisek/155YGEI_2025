from lines_to_graph2 import *
from pathlib import Path
from r2g_ver_1_2 import get_graphs, ADRESAR
import matplotlib.pyplot as plt
from kruskal import kruskal_mst
from prime import prime_mst
import heapq
from math import inf
from Bellman_ford import Bellman_ford

# def dijkstra(G, start):
#     '''
#         Dijkstra's algorithm
#         G - graph on which we are operating
#         start - the starting node
        
#         returns a list of predecessors
#     '''
    
#     # Init variables
#     n = len(G)
#     pred = [None] * (n + 1)
#     dist = [inf] * (n + 1)
#     prioQ = PriorityQueue()
    
#     # Put the start into queue and update the distance
#     prioQ.put((0, start))
#     dist[start] = 0
    
#     # While we have nodes
#     while not prioQ.empty():
#         # Get the top prio node
#         du, u = prioQ.get()
        
#         # For all it's neighbours
#         for v, wuv in G[u].items():
#             # Relax the node if there's a shorter path
#             if dist[v] > dist[u] + wuv:
#                 dist[v] = dist[u] + wuv
#                 pred[v] = u
#                 prioQ.put((dist[v], v))
                
#     return pred

def dijkstra_heap(G, start):
    n = len(G)
    dist = [inf] * n
    pred = [-1] * n

    dist[start] = 0
    heap = [(0, start)]

    while heap:
        du, u = heapq.heappop(heap)

        if du != dist[u]:   
            continue

        for v, w in G[u].items():    
            nd = du + w
            if nd < dist[v]:
                dist[v] = nd
                pred[v] = u
                heapq.heappush(heap, (nd, v))

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
    if not paths == []:
        plt.legend()
    plt.show()

def djikstra_all_pairs(G):  
    results={}

    for i, graph in enumerate(G):
        results[i]={}
        nodes=list(graph.keys())
        for s_node in nodes:
            print(s_node)
            pred=dijkstra_heap(graph, s_node)
            results[i][s_node]=pred
            # if s_node==100:
            #     return results        
    return results

def djikstra_single_pair(G, start, end):
    predecessors = [[]] * len(G)
    paths = [[]] * len(G)
    
    for i in range(len(G)):
        predecessors[i] = dijkstra_heap(G[i], start)
        
        paths[i] = rectPath(predecessors[i], start, end)
    return paths

def main():
    s_node=1
    e_node=50

    # Get graphs from OSM
    get_graphs()

    # Translate all edge data to graphs
    data_dir = Path(__file__).parent / ADRESAR
    graphs = []

    # Bonus: prepare arrays for minimum spanning trees
    kruskal_msts = []
    prime_msts = []

    for file in data_dir.iterdir():
        if file.name == "nodes.txt":
            continue
        
        new_graph = convertToGraphWithIDs(file)
        graphs.append(new_graph)

        # Make minimal spanning trees
        kruskal_msts.append(kruskal_mst(file))
        prime_msts.append(prime_mst(new_graph, 0))

        
    paths=djikstra_single_pair(graphs,s_node,e_node)

    neg_path=Bellman_ford(graphs[0],s_node)

    #=======================================================
    # Function that returns all predecessors for all nodes
    # However, the function runs for several tens of minutes up to an hour
    # Therefore, we only include the call but do not recommend running it

    # pred_all=djikstra_all_pairs(graphs)
    # a=rectPath(pred_all[0][1],1,50)
    # print(a)
    #=======================================================
    
    # Show the best paths for all weights
    # show_graph(graphs[0], paths)

    # Show the MSTs
    # for mst in kruskal_msts:
    #    show_graph(mst)

    # for mst in prime_msts:
    #     show_graph(mst)

if __name__ == "__main__":
    main()
