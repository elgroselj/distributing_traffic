from itertools import compress
import numpy as np
import networkx as nx
from scipy.optimize import linprog
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
        for counter, path in enumerate(X):
            col = edges_to_binary_vector(nodes_to_edges_path(path), list(edges))
            Q = col if Q is None else np.column_stack([Q,col])
            if counter == st_alternativ-1:
                break
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

def poklici_linprog(ZK,g,edges_mode = False):
    
    Q,B,t,c,a = sestavi_QBtca(ZK, g)

    f = np.dot(t, Q)

    # x = linprog(f=f,A=Q,b=c,Aeq=selector,beq=unit,lb=0,ub=1)
    res = linprog(f, A_ub=Q, b_ub=c, A_eq=B, b_eq=a, integrality=1)

    #print(res.fun, res.x, res.message)
    print(res.x, res.message)
    Q_used = Q[:,res.x == 1] # binarni vektorji uporabljenih poti
    paths = []
    for i in range(Q_used.shape[1]):
        if edges_mode:
            paths.append(binary_vector_to_edges(Q_used[:,i],g.edges()))
        else:
            paths.append(nodes_to_edges_path( binary_vector_to_edges(Q_used[:,i],g.edges()), inverse=True))
    return paths