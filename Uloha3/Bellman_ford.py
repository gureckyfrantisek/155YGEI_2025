from math import inf

def Bellman_ford(G, source):
    n=len(G)
    dist=[inf]*(n+1)

    for i in range(n+1):
        for edge in G:
            u,v,w=edge
            if dist[u]!=inf and dist[v]>dist[u]+w:
                
                if i==n:
                    return [-1]
                dist[v]=dist[u]+w
    return dist

