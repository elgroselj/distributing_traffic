from itertools import compress
import random
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from scipy import sparse
from collections import Counter, OrderedDict
COLORS="brgymc"

import node
import importlib
importlib.reload(node)


def plot_multigraph(graph, with_labels=True,font_size=5,figure_size=(20,20),with_arrows=True):
    plt.figure(figsize=figure_size)
    G = nx.MultiDiGraph(graph)
    try:
        pos = {n:(G.nodes[n]["x"],G.nodes[n]["y"]) for n in G.nodes()}
    except:
        
        try:
            pos = nx.planar_layout(G)
        except:
            pos = nx.circular_layout(G)
            # pos = nx.spring_layout(G)
            
    nx.draw_networkx_nodes(G, pos, node_size=0, alpha=1)
    nx.draw_networkx_labels(G, pos, font_size=font_size)
    ax = plt.gca()
    for e in G.edges:
        edge_color = G.edges[e]["color"] if "color" in G.edges[e] else "k"
        ax.annotate("",
                    xy=pos[e[1]], xycoords='data',
                    xytext=pos[e[0]], textcoords='data',
                    arrowprops=dict(arrowstyle="->" if with_arrows == True else "-", color=edge_color,
                                    shrinkA=5, shrinkB=5,
                                    patchA=None, patchB=None,
                                    connectionstyle="arc3,rad=rrr".replace('rrr', str(0.1 * (e[2] + 1))),
                                    linewidth= 7 if "over_cap" in G.edges[e] and G.edges[e]["over_cap"] else 1
                                    ),
                    )
    if with_labels:
        digraph = nx.DiGraph(G)
        nx.draw_networkx_edge_labels(digraph, pos,
                                     edge_labels={e: "{},{}".format(int(np.round(digraph.edges[e]["c"])), int(np.floor(digraph.edges[e]["cap"])))
                                                  for e in list(digraph.edges)},
                                     font_size=font_size,
                                     label_pos=0.3)

    plt.axis('off')
    plt.show()

def plot_solution_graph(graph,X,with_labels=True,font_size=5,figure_size=(20,20),with_arrows=True,over_cap_edges=[]):
    nx.set_edge_attributes(graph,{e:(e in over_cap_edges) for e in list(graph.edges())},"over_cap")
    multi = nx.MultiDiGraph(graph)
    print("k\tCOLOR")

    
    for k in range(X.shape[1]):
        path = [e for i,e in enumerate(graph.edges()) if X[i,k] != 0]
        multi.add_edges_from(path, color=COLORS[k%len(COLORS)])
        print(k,"\t",COLORS[k%len(COLORS)])
    # print(list(multi.edges()))
    plot_multigraph(multi,with_labels=with_labels,font_size=font_size,figure_size=figure_size,with_arrows=with_arrows)
    
# def plot_solution_graph_from_dict(graph,X_dict,with_labels=True,font_size=5,figure_size=(20,20)):
#     multi = nx.MultiDiGraph(graph)
#     for ie, k in X_dict:
#         multi.add_edges_from([list(graph.edges())[ie]], color=COLORS[k], label = X_dict[(ie,k)])
#     # print(list(multi.edges()))
#     plot_multigraph(multi,with_labels=with_labels,font_size=font_size,figure_size=figure_size)


##################################################################3

def fill_maxspeed(g):
    filled_neto = 0
    filled_bruto = 0
    for ke in g.edges():
        e = g.edges()[ke]
        if "maxspeed" in e.keys():
            neki = e["maxspeed"]
            if not isinstance(neki,int):
                e["maxspeed"] = int(neki) if isinstance(
                    neki, str) else int(list(neki)[0])
            else:
                print("Maxspeed already fixed.")
                return
        else:
            success = False
            for ke2 in g.edges():  # if nan set to neighbour
                if (ke[0] in ke2 or ke[1] in ke2) and ("maxspeed" in g.edges()[ke2].keys()):
                    neki = g.edges()[ke2]["maxspeed"]
                    if not isinstance(neki,int):
                        neki = int(neki) if isinstance(neki, str) else int(list(neki)[0])
                    if neki % 10 != 1:
                        e["maxspeed"] = neki
                        filled_bruto +=1
                        success = True
                        break
            if not success:
                
                def type_maxspeed(x):
                    type_speed = {"motorway":130,"trunk":110,"primary":90,"secondary":70,"tertiary":50,"unclassified":40,"residential":30}
                    key = x if isinstance(x,str) else x[0]
                    key = key.split("_")[0]
                    if key in type_speed:
                        return type_speed[key]+1
                    else:
                        return 50

                e["maxspeed"] = type_maxspeed(e["highway"])
                filled_bruto += 1
                filled_neto += 1
            
            

            print(e["highway"], end="")
            print(" is set to " + str(e["maxspeed"]))
            
    nx.set_edge_attributes(g,{e:max(10,round(g.edges()[e]["maxspeed"],-1)) for e in g.edges()},"maxspeed")
    return filled_neto, filled_bruto
            

#######################################################################################3

def nodes_to_edges_path(inp, inverse = False):
    # TODO to je dela za Q in le po cudezu za X, kaj ce mi ne pokaze vseh poti
    if not inverse:
        nl = inp
        el = []
        for i in range(1,len(nl)):
            el.append((nl[i-1], nl[i]))
        return el
    else:
        el = np.array(inp)
        #print(el)
        
        prvi_mask = ~np.isin(el[:,0], el[:,1]) 
        node = el[:,0][prvi_mask]
        
        
        #print(node)
        nl = [int(node)]
        for i in range(len(el)):
            for e in el:
                if node == e[0]:
                    node = e[1]
                    nl.append(node)
                    break
        
        return nl

def edges_to_binary_vector(el, edges):
    ev = np.zeros(len(edges))
    
    for e in el:
        ei = edges.index(e)
        if ei != -1:
            ev[ei] = 1
        
    return(ev)


def binary_vector_to_edges(ev, edges):
    return list(compress(edges, ev))
    
def columns_to_paths(g,X,edges_mode=False):
    # TODO a splitam dva avta al ne?
    paths = []
    for i in range(X.shape[1]):
        if edges_mode:
            paths.append(binary_vector_to_edges(X[:,i],g.edges()))
        else:
            paths.append(nodes_to_edges_path(binary_vector_to_edges(X[:,i],g.edges()), inverse=True))
    return paths

def split_column(col, edges):
    el = binary_vector_to_edges(col, edges)
    # TODO

def X_to_Q(g, X):
    pass # TODO
    Q = None
    for i in range(X.shape()[1]):
        new_cols = split_column(X[:,i])
        Q = new_cols if Q is None else np.column_stack([Q,new_cols])

###########################################################################################

def dfs_edges_random(graph,source,depth_limit,skip_nodes):    
    nodes = OrderedDict()
    node = source
    prev = None
    depth = 0
    while True:
        if prev is not None:
            yield (prev, node)
            
        if node not in nodes:
            successors = list(graph.successors(node))
            nodes[node] =   {
                                "successors": random.sample(successors, len(successors)),
                                "successor_idx":len(successors)-1,
                                "prev": prev,
                                "depth": depth
                            }
        prev = node
        
        b = False
        successor_idx = nodes[prev]["successor_idx"]
        if successor_idx >= 0 and depth+1 < depth_limit:
            while successor_idx >= 0:
                node = nodes[prev]["successors"][successor_idx]
                nodes[prev]["successor_idx"] -= 1
                successor_idx = nodes[prev]["successor_idx"]
                if node not in nodes and node not in skip_nodes:
                    depth += 1
                    b = True
                    break
        
        if b == False: # gremo nazaj
            if nodes[prev]["prev"] is None:
                print("prev was None")
                break
            node = nodes[prev]["prev"]
            prev = nodes[node]["prev"]
            depth -= 1
################################################333
def assemble(ks,graph,demands):
    remaining_cap = {e: graph.edges()[e]["cap"] for e in graph.edges()}
    nx.set_edge_attributes(graph, remaining_cap, "remaining_cap")
    
    alt_c = {e: graph.edges()[e]["c"] for e in graph.edges()}
    for e in graph.edges():
        if graph.edges()[e]["remaining_cap"] == 0:
            alt_c[e] = None # nasičene / s kapaciteto 0 naj bodo neskončno drage
    nx.set_edge_attributes(graph, alt_c, "alt_c")
    
    X_dict = {}
    for k in ks:              
        
        O,D,_ = demands[k]
        # time_before = time.perf_counter()
        try:
            path_of_nodes = nx.dijkstra_path(graph,O,D,"alt_c")
        except:
            print("deadend")
            return None
        # time_after = time.perf_counter()
        # Problem.print(("dijkstra time:",time_after-time_before))
        path = nodes_to_edges_path(path_of_nodes)
        
        # length = sum([graph.edges()[e]["alt_c"] for e in path])
        # if length == np.inf:
        #     # Problem.print("dead end: infinite edge was used")
        #     return None
            
        etoei = lambda e : list(graph.edges()).index(e)
        path_dict = {(etoei(e),k):1 for e in path}
        X_dict = Counter(X_dict) + Counter(path_dict)
        
        for e in path:
            graph.edges()[e]["remaining_cap"] -= 1
            if graph.edges()[e]["remaining_cap"] == 0:
                graph.edges()[e]["alt_c"] = None
    
   
    
    
    
    # Problem.print(" heuristic success")
    t = len(demands)
    m = graph.number_of_edges()
    X = sparse.dok_matrix((m,t))
    for a,k in X_dict:
        X[a,k] = X_dict[(a,k)]
    X = X.tolil()
    
    
    cap = np.array([graph.edges()[e]["cap"] for e in list(graph.edges())])
    s = cap - X.sum(axis=1).T
    if not np.all(s >= 0):
        print("false_feasible")
        return None
    return X

################################################

def latex_cell_pairs(pairs,with_field_names=True):
    repr_ = lambda x: "None" if x is None else repr(x)
    stri = r"\begin{tabular}{ll}"
    for i, (key, val) in enumerate(pairs):
        if i != 0: stri += r"\\"
        name = key if with_field_names else ""
        colon = ":" if len(name) >0 else ""
        
        stri += "{}{}         &  {}".format(name,colon,val if isinstance(val,str) else repr_(val))
            
            
    stri += "\end{tabular}"
    return stri

def latex_cell(datas):
    repr_ = lambda x: "None" if x is None else repr(x)
    stri = r"\begin{tabular}{l}"
    for i, data in enumerate(datas):
        if i != 0: stri += r"\\"
        
        stri += "{}".format(data if isinstance(data,str) else repr_(data))
            
            
    stri += "\end{tabular}"
    return stri

# def latex_row(graph,demands,):