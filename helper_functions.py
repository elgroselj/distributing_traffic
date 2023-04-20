from itertools import compress
import numpy as np
import networkx as nx
from scipy.optimize import linprog
import random
COLORS = "brgpy"


def nodes_to_edges_path(inp, inverse = False):
    if not inverse:
        nl = inp
        el = []
        for i in range(1,len(nl)):
            el.append((nl[i-1], nl[i]))
        return el
    else:
        el = np.array(inp)
        
        prvi_mask = ~np.isin(el[:,0], el[:,1]) 
        node = el[:,0][prvi_mask]
        
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

def sestavi_QBtca(ZK, G, st_alternativ = 4):
    edges = G.edges()
    edges_data = G.edges(data=True)
    Q = None
    B = None
    
    for i in range(len(ZK)):
        z,k,_ = ZK[i]
        X = nx.shortest_simple_paths(G, z, k)
        try:
            for counter, path in enumerate(X):
                col = edges_to_binary_vector(nodes_to_edges_path(path), list(edges))
                Q = col if Q is None else np.column_stack([Q,col])
                if counter == st_alternativ-1:
                    break
        except Exception as e:
            print(e)
        else:
            b = np.repeat(np.eye(len(ZK),1,-i), counter+1, axis=1)
            B = b if B is None else np.column_stack([B,b])
    
    t = [d["t"] for _,_,d in edges_data]
    c = [d["c"] for _,_,d in edges_data]
    
    a = [ai for _,_,ai in ZK]
    return(Q,B,t,c,a)

def my_nx_draw(G,paths,with_labels = False,with_nodes = True):
    pos = nx.spring_layout(G)
    if with_nodes:
        nx.draw_networkx_nodes(G,pos, node_size = 100)
    if with_labels:
        nx.draw_networkx_labels(G, pos)
    nx.draw_networkx_edges(G, pos, edge_color='k', width=0.5)
    for ci, p in enumerate(paths):
        nx.draw_networkx_edges(G, pos, edgelist=p, edge_color=COLORS[ci], width=2) # highlight elist

def poklici_linprog(ZK,g,edges_mode = False,st_alternativ=5):
    
    Q,B,t,c,a = sestavi_QBtca(ZK, g, st_alternativ= st_alternativ)
    f = np.dot(t, Q)

    # x = linprog(f=f,A=Q,b=c,Aeq=selector,beq=unit,lb=0,ub=1)
    res = linprog(f, A_ub=Q, b_ub=c, A_eq=B, b_eq=a, integrality=1)

    #print(res.fun, res.x, res.message)
    print(res.fun, res.x, res.message)
    Q_used = Q[:,res.x > 0] # binarni vektorji uporabljenih poti
    paths = []
    for i in range(Q_used.shape[1]):
        if edges_mode:
            paths.append(binary_vector_to_edges(Q_used[:,i],g.edges()))
        else:
            paths.append(nodes_to_edges_path( binary_vector_to_edges(Q_used[:,i],g.edges()), inverse=True))
    return paths


def get_random_node(G):
    nodes = list(G.nodes())
    i = random.randint(0, len(nodes)-1)
    return nodes[i]


def get_random_ZK(G,num_ZK=2,max_a=2):
    ZK = []
    for i in range(num_ZK):
        ZK.append((get_random_node(G),
                   get_random_node(G),
                   random.randint(1,max_a)))
    return ZK


def nastavi_ct(g, c_mode):
    if c_mode == "maxspeed":
        cs = {}
        for ke in g.edges():
            e = g.edges()[ke]
            if "maxspeed" in e.keys():
                neki = g.edges()[ke]["maxspeed"]
                cs[ke] = int(neki) if isinstance(
                    neki, str) else int(list(neki)[0])
            else:
                for ke2 in g.edges():  # if nan set to neighbour
                    if ke[0] in ke2 or ke[1] in ke2:
                        if "maxspeed" in g.edges()[ke2].keys():
                            neki = g.edges()[ke2]["maxspeed"]
                            cs[ke] = int(neki) if isinstance(
                                neki, str) else int(list(neki)[0])
                            g.edges()[ke]["maxspeed"] = str(cs[ke])
                            break
                if ke not in cs.keys():
                    cs[ke] = 999  # if stil set to 999

                print(e["highway"], end="")
                print(" is set to " + str(cs[ke]))
        nx.set_edge_attributes(g, cs, name='c')
        nx.set_edge_attributes(
            g, {e: g.edges()[e]["length"]/g.edges()[e]["c"] for e in g.edges()}, name='t')
    elif isinstance(c_mode, int):
        nx.set_edge_attributes(g, {e: c_mode for e in g.edges()}, name='c')
        nx.set_edge_attributes(
            g, {e: g.edges()[e]["length"] for e in g.edges()}, name='t')
