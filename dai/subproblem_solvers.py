from collections import Counter
import networkx as nx
import numpy as np
from scipy import sparse

def dijkstra1(vp, graph, demands, alpha):
    X_dict = {}
    status = "feasible"
    
    _,m = vp["B"].shape 
    _,t = vp["H"].shape                
    new_c = vp["c"].value + vp["lam"].value
    alt_c = {e: new_c[ei] for ei, e in enumerate(graph.edges())} 
    nx.set_edge_attributes(graph, alt_c, "alt_c")
    
    for k in range(t):
        Ok, Dk, num_k = demands[k]
        
        On = list(graph.nodes())[Ok]
        Dn = list(graph.nodes())[Dk]
        
        try:
            path_of_nodes = nx.dijkstra_path(graph,On,Dn,"alt_c")
        except:
            print("path does not exist")
            status = "infeasible"
            break
        
        def nodes_to_edges_path(inp):
            nl = inp
            el = []
            for i in range(1,len(nl)):
                el.append((nl[i-1], nl[i]))
            return el
        path = nodes_to_edges_path(path_of_nodes)
        # length = sum([graph.edges()[e]["alt_c"] for e in path])
        # print(length)
        etoei = lambda e : list(graph.edges()).index(e)
        path_dict = {(etoei(e),k):num_k for e in path}
        X_dict = Counter(X_dict) + Counter(path_dict)
        
    X = sparse.dok_matrix((m,t))
    for a,k in X_dict:
        X[a,k] = X_dict[(a,k)]
    X = np.array(X.todense())
    
    zLD = new_c.T @ np.sum(X,axis=1) - vp["lam"].value.T @ vp["cap"].value
    zLD = int(zLD[0])
    s = vp["cap"].value - np.sum(X,axis=1)
    
    ll = vp["lam"].value - s * alpha
    ll[ll < 0] = 0 # lambda ne mora biti negativna
    vp["lam"].value = ll
    
    return (X, status, zLD, s)

def dijkstra2(vp, graph, demands, U, L, alpha):
    X_dict = {}
    status = "feasible"
    
    _,m = vp["B"].shape 
    _,t = vp["H"].shape
    

    
    for k in range(t):
        new_c = vp["c"].value + vp["lam"].value + vp["eta"].value[:,k] - vp["xi"].value[:,k]
        new_c[new_c < 0] = 0
        alt_c = {e: new_c[ei] for ei, e in enumerate(graph.edges())} 
        nx.set_edge_attributes(graph, alt_c, "alt_c")
        
        
        Ok, Dk, num_k = demands[k]
        
        On = list(graph.nodes())[Ok]
        Dn = list(graph.nodes())[Dk]
        
        try:
            path_of_nodes = nx.dijkstra_path(graph,On,Dn,"alt_c")
        except:
            print("path does not exist")
            status = "infeasible"
            break
        
        def nodes_to_edges_path(inp):
            nl = inp
            el = []
            for i in range(1,len(nl)):
                el.append((nl[i-1], nl[i]))
            return el
        path = nodes_to_edges_path(path_of_nodes)
        # length = sum([graph.edges()[e]["alt_c"] for e in path])
        # print(length)
        etoei = lambda e : list(graph.edges()).index(e)
        path_dict = {(etoei(e),k):num_k for e in path}
        X_dict = Counter(X_dict) + Counter(path_dict)
        
    X = sparse.dok_matrix((m,t))
    for a,k in X_dict:
        X[a,k] = X_dict[(a,k)]
    X = np.array(X.todense())
    
    zLD = new_c.T @ np.sum(X, axis=1) - vp["lam"].value.T @ vp["cap"].value - np.sum(
                            vp["eta"].value.T @ U) - np.sum(vp["xi"].value.T @ L)
    # zLD = int(zLD[0])
    s_lam = vp["cap"].value - np.sum(X, axis=1)
    s_eta = U - X
    s_xi = X - L
    
    ll = vp["lam"].value - s_lam * alpha
    ll[ll < 0] = 0 # lambda ne mora biti negativna
    vp["lam"].value = ll
    
    ee = vp["eta"].value - s_eta * alpha
    ee[ee < 0] = 0 # lambda ne mora biti negativna
    vp["eta"].value = ee
    
    xx =  vp["xi"].value - s_xi * alpha
    xx[xx < 0] = 0 # lambda ne mora biti negativna
    vp["xi"].value = xx
    
    s = np.concatenate((s_lam,s_eta.flatten(),s_xi.flatten()))
    return (X, status, zLD, s)

def cvxpy_linprog_LR(problem, vp,alpha):
    problem.solve()
    
    X = vp["X"].value
    status = problem.status
    if status == "infeasible":
        return (None,status,None,None)
    zLD = problem.value
    s = vp["cap"].value - np.sum(X,axis=1)
    
    ll = vp["lam"].value - s * alpha
    ll[ll < 0] = 0 # lambda ne mora biti negativna
    vp["lam"].value = ll
    
    return (X, status, zLD, s)