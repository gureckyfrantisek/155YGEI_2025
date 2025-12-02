# function to download and process the graph from OSM

import os
import osmnx as ox
import math
import re  

MISTO = 'Moravskoslezský kraj, Czechia' 
ADRESAR = 'data_grafy' 

# File names
FILE_NODES = os.path.join(ADRESAR, 'nodes.txt')
FILE_DIST = os.path.join(ADRESAR, 'graph_dist.txt')
FILE_TIME1 = os.path.join(ADRESAR, 'graph_time_basic.txt')
FILE_TIME2 = os.path.join(ADRESAR, 'graph_time_tortuous.txt')
FILE_TORT = os.path.join(ADRESAR, 'graph_tortuosity.txt')

def get_graphs():
    ALL_FILES = [FILE_NODES, FILE_DIST, FILE_TIME1, FILE_TIME2, FILE_TORT]

    if not os.path.exists(ADRESAR):
        os.makedirs(ADRESAR)

    # CHECK IF FILES EXIST
    files_exist = all(os.path.exists(f) for f in ALL_FILES)

    if files_exist:
        print(f"Files already exist in folder '{ADRESAR}'.")
        return

    ox.settings.timeout = 2000
    # Download graph
    G_raw = ox.graph_from_place(MISTO, network_type='drive')
    # Project and convert to undirected graph
    G_proj = ox.project_graph(G_raw)
    G = ox.convert.to_undirected(G_proj)

    # Compute weights
    edges_dist = []
    edges_time1 = []
    edges_time2 = []
    edges_tort = []

    # Dictionary to convert textual speeds to numbers
    SPEED_DICT = {
        'CZ:urban': 50,
        'CZ:rural': 90,
        'CZ:trunk': 110,
        'CZ:motorway': 130,
        'urban': 50,
        'rural': 90,
        'none': 130,
        'walk': 5,
        'living_street': 20
    }

    # Helper dict to map OSM node ids to consecutive ids (0..n)
    # Format OSM_id: new_id
    ordered_nodes = {}
    nodes_counter = 0

    for u, v, data in G.edges(data=True):
        # Distance
        length_m = data['length']
        
        # Speed
        max_speed_raw = data.get('maxspeed', '50')
        
        # If it's a list, take the first value
        if isinstance(max_speed_raw, list):
            max_speed_raw = max_speed_raw[0]
            
        max_speed_str = str(max_speed_raw).strip()
        
        # Try to find in the translation dictionary
        if max_speed_str in SPEED_DICT:
            final_speed = SPEED_DICT[max_speed_str]
        else:
            # If not in dict, try to extract a number from the string
            # Use regex to find the first number in the text (e.g., "50 mph" -> 50)
            found_numbers = re.findall(r'\d+', max_speed_str)
            if found_numbers:
                try:
                    final_speed = int(found_numbers[0])
                except ValueError:
                    final_speed = 50 
            else:
                final_speed = 50 

        # Handle zero or negative speed
        if final_speed <= 0:
            final_speed = 50

        speed_mps = final_speed / 3.6

        # Tortuosity
        x_u, y_u = G.nodes[u]['x'], G.nodes[u]['y']
        x_v, y_v = G.nodes[v]['x'], G.nodes[v]['y']
        euclid_dist = math.sqrt((x_u - x_v)**2 + (y_u - y_v)**2)
        
        if euclid_dist < 0.1:
            tortuosity = 1.0
        else:
            tortuosity = length_m / euclid_dist
        
        # Times
        time_basic = length_m / speed_mps
        
        effective_speed = speed_mps / tortuosity
        if effective_speed > 0:
            time_tortuous = length_m / effective_speed 
        else:
            time_tortuous = 999999

        # Save in format: u;v;weight

        # Check if node was already used, if not assign a new id
        if u not in ordered_nodes:
            ordered_nodes[u] = nodes_counter
            nodes_counter += 1
        
        if v not in ordered_nodes:
            ordered_nodes[v] = nodes_counter
            nodes_counter += 1
        
        u = ordered_nodes[u]
        v = ordered_nodes[v]

        # Edge entries
        edges_dist.append(f"{u} {v} {length_m:.2f}")
        edges_time1.append(f"{u} {v} {time_basic:.2f}")
        edges_time2.append(f"{u} {v} {time_tortuous:.2f}")
        edges_tort.append(f"{u} {v} {tortuosity:.4f}")

    # Saving to files

    with open(FILE_NODES, 'w') as f:
        for node, data in G.nodes(data=True):
            f.write(f"{ordered_nodes[node]} {data['x']:.2f} {data['y']:.2f}\n")

    def save_edges(filename, lines):
        with open(filename, 'w') as f:
            for line in lines:
                f.write(line + "\n")

    save_edges(FILE_DIST, edges_dist)
    save_edges(FILE_TIME1, edges_time1)
    save_edges(FILE_TIME2, edges_time2)
    save_edges(FILE_TORT, edges_tort)