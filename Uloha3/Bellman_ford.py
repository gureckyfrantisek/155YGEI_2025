from math import inf

def Bellman_ford(G, start):
    nodes=list(G.keys())
    n=len(nodes)
    dist={node:inf for node in nodes}
    dist[start]=0
    for _ in range(n-1):
        updated=False
        for u in G:
            if dist[u]==inf:
                continue
            for v,w in G[u].items():
                if dist[v]>dist[u]+w:
                    dist[v]=dist[u]+w
                    updated=True
        if not updated:
            break
    return dist

