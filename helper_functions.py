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

def edge_to_incidence_vector(e, nodes):
    nv = np.zeros(len(nodes))
    
    ni = nodes.index(e[0])
    nv[ni] = 1
    ni = nodes.index(e[1])
    nv[ni] = -1
        
    return(nv)

def binary_vector_to_edges(ev, edges):
    return list(compress(edges, ev))

def sestavi_QBtca(ZK, G, st_alternativ = 4):
    edges = G.edges()
    edges_data = G.edges(data=True)
    Q = None
    B = None
    
    for i in range(len(ZK)):
        z,k,_ = ZK[i]
        X = nx.shortest_simple_paths(G, z, k, weight="t")
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

def columns_to_paths(g,X,edges_mode=False):
    paths = []
    for i in range(X.shape[1]):
        if edges_mode:
            paths.append(binary_vector_to_edges(X[:,i],g.edges()))
        else:
            paths.append(nodes_to_edges_path( binary_vector_to_edges(X[:,i],g.edges()), inverse=True))
    return paths

def poklici_linprog(ZK,g,edges_mode = False,st_alternativ=5,integrality=1):
    
    Q,B,t,c,a = sestavi_QBtca(ZK, g, st_alternativ= st_alternativ)
    f = np.dot(t, Q)

    # x = linprog(f=f,A=Q,b=c,Aeq=selector,beq=unit,lb=0,ub=1)
    res = linprog(f, A_ub=Q, b_ub=c, A_eq=B, b_eq=a, integrality=1)

    #print(res.fun, res.x, res.message)
    print(res.fun, res.x, res.message)
    Q_used = Q[:,res.x > 0] # binarni vektorji uporabljenih poti
    paths = columns_to_paths(g,Q_used,edges_mode=edges_mode)
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

def fill_maxspeed(g):
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
                    if neki != 999:
                        e["maxspeed"] = neki
                        success = True
                        break
            if not success:
                e["maxspeed"] = 999

            print(e["highway"], end="")
            print(" is set to " + str(e["maxspeed"]))

def nastavi_c(g, c_mode=None, c=None):
    if c is not None:
        cs = {e : c[i] for i, e in enumerate(g.edges())}
    elif c_mode == "maxspeed":
        fill_maxspeed(g)
        maxspeed_to_capacity = lambda maxspeed: int(maxspeed/10)
        cs = {e : maxspeed_to_capacity(g.edges()[e]["maxspeed"]) for e in g.edges()}
    elif isinstance(c_mode, int):
        cs = {e: c_mode for e in g.edges()}
    else:
        raise Exception("Invalid mode or c missing.")
    nx.set_edge_attributes(g, cs, name='c')

def nastavi_t(g, t_mode=None, t=None):
    if t is not None:
        ts = {e : t[i] for i, e in enumerate(g.edges())}
    elif t_mode == "length":
        ts = {e: g.edges()[e]["length"] for e in g.edges()}
    elif t_mode == "time":
        fill_maxspeed(g)
        ts = {e: g.edges()[e]["length"]/g.edges()[e]["maxspeed"] for e in g.edges()}
    elif isinstance(t_mode, int):
        ts = {e: t_mode for e in g.edges()}
    else:
        raise Exception("Invalid mode or t missing.")
    nx.set_edge_attributes(g, ts, name='t')
        
def nastavi_ct(g, c_mode, c = None, t = None):
    if c_mode == "maxspeed":
        fill_maxspeed(g)
        maxspeed_to_capacity = lambda maxspeed: int(maxspeed/10)
        cs = {e : maxspeed_to_capacity(g.edges()[e]["maxspeed"]) for e in g.edges()}
            
        nx.set_edge_attributes(g, cs, name='c')
        nx.set_edge_attributes(
            g, {e: g.edges()[e]["length"]/g.edges()[e]["maxspeed"] for e in g.edges()}, name='t')
    elif isinstance(c_mode, int):
        nx.set_edge_attributes(g, {e: c_mode for e in g.edges()}, name='c')
        nx.set_edge_attributes(
            g, {e: g.edges()[e]["length"] for e in g.edges()}, name='t')

from scipy.sparse import coo_array
def sparse_incidence_matrix(nodes, edges, factor = None):
    nodes = list(nodes)
    edges = list(edges)
    row_nodes = []
    col_edges = []
    data = []
    for ei, e in enumerate(edges):
        if factor is not None:
            k = factor[ei]
        else:
            k = 1
        ni = nodes.index(e[0])
        row_nodes.append(ni)
        col_edges.append(ei)
        data.append(1*k)
        ni = nodes.index(e[1])
        row_nodes.append(ni)
        col_edges.append(ei)
        data.append(-1*k)

    coo = coo_array((data, (row_nodes, col_edges)), shape=(len(nodes), len(edges)))
    return coo

def getB(num_edges, num_ZK, factors = None):
    B = None
    row_edge = []
    col_ZK = []
    data = []
    q = 0
    for i in range(num_edges):
        for j in range(num_ZK):
            row_edge.append(i)
            col_ZK.append(q)
            if factors is None:
                data.append(1)
            else:
                data.append(factors[i])
            q += 1
            
            
        # b = np.repeat(np.eye(num_edges,1,-i), num_ZK, axis=1)
        # B = b if B is None else np.column_stack([B,b])
    B = coo_array((data, (row_edge, col_ZK)), shape=(num_edges, num_ZK*num_edges))
    return B

def MtoM2(M,num_ZK):
    num_nodes, num_edges = M.shape
    # TODO
    rows = []
    cols = []
    data = []
    
    shift = 0
    for k in range(num_ZK):
        
        for n in range(num_nodes):
            for e in range(num_edges):
                rows.append(k*num_nodes+n)
                cols.append(e*num_ZK + shift)
                data.append(M[n,e])
        shift += 1
    coo = coo_array((data, (rows, cols)), shape=(num_ZK*num_nodes, num_ZK*num_edges))
    return coo

def nastavi_fbcmm(G,ZK,c,t):
    f = np.repeat(t,len(ZK)) #min f*x
    
    B = getB(G.number_of_edges(),len(ZK)) # B*x <= c
    
    M = sparse_incidence_matrix(G.nodes(),G.edges()) # (M*X=M_ZK)
    
    a = [a for _,_,a in ZK]
    M_ZK = sparse_incidence_matrix(G.nodes(),[(z,k) for z, k, _ in ZK],factor=a)
    
    M2 = MtoM2(M.todense(), len(ZK)) # M2 * x = m_ZK
    m_ZK = M_ZK.todense().flatten('F')
    
    return (f,B,c,M2,m_ZK)