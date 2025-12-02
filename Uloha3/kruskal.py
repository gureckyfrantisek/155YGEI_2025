# Read the file again
def read_edges_with_sort(file):
    # If we've read the file before, set the pointer to the begining
    if hasattr(file, "read"):
        f = file
        try:
            f.seek(0)
        except Exception:
            pass
    else:
        f = open(file, "r")

    edges = []
    for line in f:
        sid, eid, w = line.split()
        edges.append((int(sid), int(eid), float(w)))

    f.close()

    # Sort the edges by weight
    edges.sort(key=lambda edge: edge[2])
    return edges

# Weighted Union and Find with path compression from lecture
def find(u, p):
    while p[u] != u:
        # If parent is not the root, move to their grandparent and save a step
        u = p[p[u]]
        u = p[u]
    return u

def union(u, v, p, r):
    root_u = find(u, p)
    root_v = find(v, p)

    if root_u != root_v:
        # If not in the same subtree
        # Connect lower ranked to higher rank
        if r[root_u] > r[root_v]:
            p[root_v] = root_u
        
        elif r[root_v] > r[root_u]:
            p[root_u] = root_v

        # If equal ranks, connect together and increment rank
        else:
            p[root_u] = root_v
            r[root_v] = r[root_v] + 1

        return True
    
    return False

def kruskal_mst(file):
    edges = read_edges_with_sort(file)

    # Collect nodes, use set for quicker search
    nodes = set()
    for u, v, _ in edges:
        nodes.add(u)
        nodes.add(v)

    parent = {n: n for n in nodes}
    rank = {n: 0 for n in nodes}

    MST = {n: {} for n in nodes}
    for u, v, w in edges:
        if union(u, v, parent, rank):
            MST[u][v] = w
            MST[v][u] = w

    return MST