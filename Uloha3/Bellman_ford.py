from math import inf

def Bellman_ford(G, start):
    nodes=list(G.keys())
    n=len(nodes)
    dist={node:inf for node in nodes}
    dist[start]=0
    pred={node:None for node in nodes}
    for i in range(n):
        change=False
        for u in G:
            if dist[u]==inf:
                continue
            for v,w in G[u].items():
                if dist[v]>dist[u]+w:
                    dist[v]=dist[u]+w
                    pred[v]=u
                    change=True
        if not change:
            break
        if i==n-1:
            raise ValueError("Graph has negative weight cycle")

    # # Check for negative weight cycles
    # for u in G:
    #     if dist[u]==inf:
    #         continue
    #     for v,w in G[u].items():
    #         if dist[v]>dist[u]+w:
    #             raise ValueError("Graph has negative weight cycle")
    return pred, dist

