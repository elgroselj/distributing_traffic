from matplotlib import pyplot as plt
import networkx as nx
import osmnx as ox
import random

import numpy as np
import helper_functions as hf
import pickle


def relabel(graph, demands):
    graph = nx.convert_node_labels_to_integers(
        graph, label_attribute="old_name")

    def f(old_name): return [n for (n, d) in graph.nodes(
        data=True) if d["old_name"] == old_name][0]
    demands = [(f(O), f(D), num_k) for (O, D, num_k) in demands]
    return graph, demands

# def smooth_out_simple_nodes(graph_,demands,min_c):
#     nodes_to_remove = [n for n in graph_.nodes if len(list(graph_.predecessors(n))) == 1 and len(list(graph_.successors(n))) == 1]

#     graph = nx.MultiDiGraph(graph_)
#     nodes_from_demands = sum([[O,D] for (O,D,_) in demands],[])
#     # For each of those nodes
#     for node in nodes_to_remove:
#         if node in nodes_from_demands: # do not remove
#             continue

#         in_node = list(graph.predecessors(node))[0]
#         out_node = list(graph.successors(node))[0]

#         if len(graph[in_node][node]) > 1 or len(graph[node][out_node]) > 1:
#             continue

#         c_sum = graph[in_node][node][0]["c"] + graph[node][out_node][0]["c"]
#         if c_sum > min_c:
#             continue

#         new_cap = min(graph[in_node][node][0]["cap"], graph[node][out_node][0]["cap"])

#         graph.add_edge(in_node, out_node, c = c_sum, cap = new_cap)
#         graph.remove_node(node)

#     graph,demands = relabel(graph,demands)

#     return (graph, demands)


def smooth_out_simple_nodes(graph_, demands, min_c):
    graph = nx.DiGraph(graph_)
    nodes_to_remove = [n for n in graph.nodes if len(
        list(graph.predecessors(n))) == 1 and len(list(graph.successors(n))) == 1]

    nodes_from_demands = sum([[O, D] for (O, D, _) in demands], [])
    # For each of those nodes
    for node in nodes_to_remove:
        if node in nodes_from_demands:  # do not remove
            continue

        in_node = list(graph.predecessors(node))[0]
        out_node = list(graph.successors(node))[0]

        if (in_node, out_node) in graph.edges():  # nočem multigrapha
            continue

        new_c = graph[in_node][node]["c"] + graph[node][out_node]["c"]
        if new_c > min_c:
            continue

        new_cap = min(graph[in_node][node]["cap"],
                      graph[node][out_node]["cap"])

        graph.add_edge(in_node, out_node, c=new_c, cap=new_cap)
        graph.remove_node(node)

    graph, demands = relabel(graph, demands)

    return (graph, demands)

    # nx.draw(G,with_labels=True)


# ta density loh spreminjam če bojo vsi enako gosti
def generate_random_graph(n_max, t, k_num_max, cap_max, density_param=0.5, c_max=100):

    # generiraš podatke

    edges = set()
    demands = []
    for k in range(t):
        O = random.randint(0, n_max-1)
        D = random.randint(0, n_max-2)
        if D >= O:
            D += 1   # O != D

        k_num = random.randint(1, k_num_max)
        demands.append((O, D, k_num))
        for ki in range(k_num):
            # (k_num_max/2) je pričakovana vrednost
            len_ub = min(n_max-1, density_param *
                         (n_max-1)*n_max/t/(k_num_max/2))
            len = random.randint(1, len_ub)
            nodes = random.sample(
                [node for node in range(n_max) if node not in [O, D]], len-1)
            edges_ = hf.nodes_to_edges_path([O]+nodes+[D])
            edges.update(edges_)

    edges = list(edges)

    graph = nx.DiGraph()
    graph.add_edges_from(edges)

    c = {e: random.randint(1, c_max) for e in list(graph.edges())}
    cap = {e: random.randint(1, cap_max) for e in list(graph.edges())}

    # konstruiraš graf iz podatkov

    nx.set_edge_attributes(graph, c, "c")
    nx.set_edge_attributes(graph, cap, "cap")

    graph, demands = relabel(graph, demands)

    print("demands: ", demands)
    # hf.plot_multigraph(graph, with_labels=True,font_size=10,figure_size=(20,20))

    return (graph, demands)

# graph, demands = generate_random_graph(n_max=10,t=2,k_num_max=2,cap_max=2,density_param=1)
# # graph = nx.DiGraph(nx.house_graph())

# print([x for x in nx.dfs_edges(graph,source=0,depth_limit=3)])
# #

# # print(demands)
# # print(nx.density(graph))


# hf.plot_multigraph(graph, with_labels=True,font_size=10,figure_size=(20,20))


def place_to_nx(place, save=False, mode="place", point=None):
    try:
        file = open('/home/lema/Documents/diplomska/dai/' +
                    place+'_nx.pkl', 'rb')
        graph = pickle.load(file)
        file.close()
        return graph
    except:
        pass

    if mode == "place":
        # graph_ox = ox.graph_from_place(place,network_type='drive',custom_filter='["highway"~"primary"]')
        graph_ox = ox.graph_from_place(place, network_type='drive')
    elif mode == "point":
        # (46.05407,14.52114)
        graph_ox = ox.graph_from_point(point, dist=10000, network_type="drive")

    graph = nx.convert_node_labels_to_integers(nx.DiGraph(graph_ox))

    hf.fill_maxspeed(graph)

    graph.edges(data=True)
    times = {e: np.round(
        graph.edges()[e]["length"]/graph.edges()[e]["maxspeed"]) for e in graph.edges()}
    capacities = {e: np.floor(
        1 + graph.edges()[e]["length"]*graph.edges()[e]["maxspeed"]/1000) for e in graph.edges()}
    nx.set_edge_attributes(graph, times, "c")
    nx.set_edge_attributes(graph, capacities, "cap")

    if save:
        file = open('/home/lema/Documents/diplomska/dai/' +
                    place+'_nx.pkl', 'wb')
        pickle.dump(graph, file)
        file.close()

    return graph

# place = 'Castenaso'
# graph = place_to_nx(place,save=True)


def point_to_closest_node(point, graph):
    f=lambda node_and_data: abs(node_and_data[1]["y"] - point[0]) + abs(node_and_data[1]["x"] - point[1])
        
    closest_node, _ = min([(node,data) for node,data in graph.nodes(data=True)],
                       key=lambda node_and_data: f(node_and_data))
    return closest_node


def random_demands(graph,num_cars):
    demands = []
    for car in range(num_cars):
        while True:
            O,D = random.sample(graph.nodes(),2)
            try:
                nx.dijkstra_path(graph,O,D)
            except:
                continue
            demands.append((O,D,1))
            break
    unique_OD = set([(O,D) for O,D,_ in demands])
    demands_ = []
    for O_,D_ in unique_OD:
        demands_.append((O_,D_,len([1 for O,D,_ in demands if (O == O_ and D == D_)])))
    return demands_
        
    
        
        
    return demands

def about_graph(graph):
    print("graph",end=": ")
    print("num_nodes: ", graph.number_of_nodes(),end=", ")
    print("density: ", nx.density(graph))
    pass
    
    
def about_demands(demands):
    print("demands: ",demands,end=": ")
    print("num_cars: ", sum([k_num for _,_,k_num in demands]),end=", ")
    print("avg_group_size: ", sum([k_num for _,_,k_num in demands])/len(demands))
    print("biggest_group_size: ", max([k_num for _,_,k_num in demands]))
