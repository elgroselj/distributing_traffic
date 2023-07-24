import networkx as nx
import osmnx as ox
import random

import numpy as np
import helper_functions as hf
import pickle

def generate_random_graph(n_max,t,k_num_max,cap_max,density_param=0.5): # ta density loh spreminjam če bojo vsi enako gosti
    
    # generiraš podatke
    
    edges = set()
    demands = []
    for k in range(t):
        O = random.randint(0,n_max-1)
        D = random.randint(0,n_max-2)
        if D >= O: D += 1   # O != D
        
        k_num = random.randint(1,k_num_max)
        demands.append((O,D,k_num))
        for ki in range(k_num):
            len_ub = min(n_max-1,density_param*(n_max-1)*n_max/t/(k_num_max/2)) # (k_num_max/2) je pričakovana vrednost
            len = random.randint(1,len_ub)
            nodes = random.sample([node for node in range(n_max) if node not in [O,D]],len-1)
            edges_ = hf.nodes_to_edges_path([O]+nodes+[D])
            edges.update(edges_)
    
    edges = list(edges)
    
    
    graph = nx.DiGraph()
    graph.add_edges_from(edges)
    
    C_MAX = 100
    c = {e : random.randint(1,C_MAX) for e in list(graph.edges())}
    cap = {e : random.randint(1,cap_max) for e in list(graph.edges())}
    
    # konstruiraš graf iz podatkov
    
 
    
    
    nx.set_edge_attributes(graph, c,"c")
    nx.set_edge_attributes(graph, cap,"cap")
    
    
    graph = nx.convert_node_labels_to_integers(graph,label_attribute="old_name")
    f = lambda old_name: [n for (n,d) in graph.nodes(data=True) if d["old_name"] == old_name][0]
    demands = [(f(O),f(D),num_k) for (O,D,num_k) in demands]
    
    print("demands: ", demands)
    hf.plot_multigraph(graph, with_labels=True,font_size=10,figure_size=(20,20))
    
    return (graph, demands)
              
    
    
# 

# print(demands)
# print(nx.density(graph))


# hf.plot_multigraph(graph, with_labels=True,font_size=10,figure_size=(20,20))

    


def place_to_nx(place,save = False):
    try:
        file = open('/home/lema/Documents/diplomska/dai/'+place+'_nx.pkl', 'rb')
        graph = pickle.load(file)
        file.close()
        return graph
    except:
        pass
    
    graph_ox = ox.graph_from_place(place,network_type='drive')

    graph = nx.convert_node_labels_to_integers(nx.DiGraph(graph_ox))

    hf.fill_maxspeed(graph)

    graph.edges(data=True)
    times = {e: np.round(graph.edges()[e]["length"]/graph.edges()[e]["maxspeed"]) for e in graph.edges()}
    capacities = {e: np.floor(1 + graph.edges()[e]["length"]*graph.edges()[e]["maxspeed"]/1000) for e in graph.edges()}
    nx.set_edge_attributes(graph, times,"c")
    nx.set_edge_attributes(graph, capacities,"cap")
    
    if save:
        file = open('/home/lema/Documents/diplomska/dai/'+place+'_nx.pkl', 'wb')
        pickle.dump(graph,file)
        file.close()
    
    return graph
    
# place = 'Castenaso'
# graph = place_to_nx(place,save=True)