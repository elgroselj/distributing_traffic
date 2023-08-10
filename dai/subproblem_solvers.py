from collections import Counter
import random
import networkx as nx
import numpy as np
from scipy import sparse
import helper_functions as hf

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
    # zLD = int(zLD[0])
    s = vp["cap"].value - np.sum(X,axis=1)
    
    ll = vp["lam"].value - s * alpha
    ll[ll < 0] = 0 # lambda ne mora biti negativna
    vp["lam"].value = ll
    
    return (X, status, zLD, s)



def dijkstra11(vp, graph, demands, alpha, without_edges):
    X_dict = {}
    status = "feasible"
    
    _,m = vp["B"].shape 
    _,t = vp["H"].shape                
    new_c = vp["c"].value + vp["lam"].value
    alt_c = {e: new_c[ei] for ei, e in enumerate(graph.edges())}
    
    for ei in without_edges:
        e = list(graph.edges())[ei]
        alt_c[e] = np.inf
        
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
        
        path = hf.nodes_to_edges_path(path_of_nodes)
        length = sum([graph.edges()[e]["alt_c"] for e in path])
        if length == np.inf:
            # tj uporabili smo prepovedano povezavo
            X = None
            status = "infeasible"
            zLD = None
            s = None
            return (X, status, zLD, s)
        # print(length)
        etoei = lambda e : list(graph.edges()).index(e)
        path_dict = {(etoei(e),k):num_k for e in path}
        X_dict = Counter(X_dict) + Counter(path_dict)
        
    X = sparse.dok_matrix((m,t))
    for a,k in X_dict:
        X[a,k] = X_dict[(a,k)]
    X = np.array(X.todense())
    
    zLD = new_c.T @ np.sum(X,axis=1) - vp["lam"].value.T @ vp["cap"].value
    # zLD = int(zLD[0])
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

def dijkstra0(vp, graph, demands, alpha):
    X_dict = {}
    status = "neki"
    
    _,m = vp["B"].shape 
    _,t = vp["H"].shape                
    new_c = vp["c"] + vp["lam"]
    alt_c = {e: new_c[0,ei] for ei, e in enumerate(graph.edges())} 
    nx.set_edge_attributes(graph, alt_c, "alt_c")
    
    for k in range(t):
        Ok, Dk, num_k = demands[k]
        
        On = list(graph.nodes())[Ok]
        Dn = list(graph.nodes())[Dk]
        
        try:
            path_of_nodes = nx.dijkstra_path(graph,On,Dn,"alt_c")
        except:
            # print("path does not exist")
            status = "infeasible"
            break
        
        path = hf.nodes_to_edges_path(path_of_nodes)
        # length = sum([graph.edges()[e]["alt_c"] for e in path])
        # print(length)
        etoei = lambda e : list(graph.edges()).index(e)
        path_dict = {(etoei(e),k):num_k for e in path}
        X_dict = Counter(X_dict) + Counter(path_dict)
        
    X = sparse.dok_matrix((m,t))
    for a,k in X_dict:
        X[a,k] = X_dict[(a,k)]
    X = X.tolil()
    
    # zLD = new_c.T @ np.sum(X,axis=1) - vp["lam"].T @ vp["cap"]
    zLD = (new_c@X.sum(axis=1))[0,0] - (vp["lam"] @ vp["cap"])[0]
    # zLD = int(zLD[0])
    s = vp["cap"] - X.sum(axis=1).T
    
    ll = vp["lam"] - s * float(alpha)
    ll = vp["lam"] - s * float(alpha)
    ll[ll < 0] = 0 # lambda ne mora biti negativna
    vp["lam"] = ll
    
    return (X, status, zLD, s)

def astar(vp, graph, demands, alpha):
    
    X_dict = {}
    status = "feasible"
    
    _,m = vp["B"].shape 
    _,t = vp["H"].shape                
    new_c = vp["c"] + vp["lam"]
    alt_c = {e: new_c[0,ei] for ei, e in enumerate(graph.edges())} 
    nx.set_edge_attributes(graph, alt_c, "alt_c")
    
    for k in range(t):
        Ok, Dk, num_k = demands[k]
        
        On = list(graph.nodes())[Ok]
        Dn = list(graph.nodes())[Dk]
        
        try:
            # path_of_nodes = nx.dijkstra_path(graph,On,Dn,"alt_c")
            dist = lambda a,b: ((graph.nodes()[a]["x"] - graph.nodes()[b]["x"])**2 + (graph.nodes()[a]["y"] - graph.nodes()[b]["y"])**2)**0.5
            path_of_nodes = nx.astar_path(graph,On,Dn,dist,"alt_c")
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
    X = X.tolil()
    
    # zLD = new_c.T @ np.sum(X,axis=1) - vp["lam"].T @ vp["cap"]
    zLD = (new_c@X.sum(axis=1))[0,0] - (vp["lam"] @ vp["cap"])[0]
    # zLD = int(zLD[0])
    s = vp["cap"] - X.sum(axis=1).T
    
    ll = vp["lam"] - s * float(alpha)
    ll = vp["lam"] - s * float(alpha)
    ll[ll < 0] = 0 # lambda ne mora biti negativna
    vp["lam"] = ll
    
    return (X, status, zLD, s)

def epaths_to_X(graph,epaths,ks,demands):
    X_dict = {}
    for ki, k in enumerate(ks):
        etoei = lambda e : list(graph.edges()).index(e)
        path_dict = {(etoei(e),k):1 for e in epaths[ki]}
        X_dict = Counter(X_dict) + Counter(path_dict)
    
    m = graph.number_of_edges()
    t = len(demands)
    X = sparse.dok_matrix((m,t))
    for a,k in X_dict:
        X[a,k] = X_dict[(a,k)]
    X = X.tolil()
    return X

def dijkstra0_with_without(vp, graph, demands, alpha, with_edges=[],without_edges=[]):
    X = None
    status = "neki"
    
    _,m = vp["B"].shape 
    _,t = vp["H"].shape                
    new_c = vp["c"] + vp["lam"].value
    try:
        alt_c = {e: (None if (e in without_edges) else new_c[ei]) for ei, e in enumerate(graph.edges())} 
    except:
        pass
    nx.set_edge_attributes(graph, alt_c, "alt_c")
    
    
    epaths = []
    for k in range(t):
        Ok, Dk, num_k = demands[k]
        
        On = list(graph.nodes())[Ok]
        Dn = list(graph.nodes())[Dk]
        
        try:
            path_of_nodes = nx.dijkstra_path(graph,On,Dn,"alt_c")
        except:
            # print("path does not exist")
            status = "infeasible"
            return None, "infeasible", None, None
        
        path = hf.nodes_to_edges_path(path_of_nodes)
        epaths += [path] * num_k
       
    
    ks = sum([[i]*data[2] for i,data in enumerate(demands)],[])
    free_ks = []
    for ki, k in enumerate(ks):
        inters = set(epaths[ki]).intersection(with_edges)
        if len(inters) == 0:
            free_ks.append((ki,k))
        else:
            # removed covered edges
            with_edges = [e for e in with_edges if e not in inters]
            
    if len(with_edges) == 0:
        pass
    
    elif len(with_edges) > len(free_ks):
        raise Exception("Too many with edges")
    
    else:
    
        for edge in with_edges:
            old_alt_c = graph.edges()[edge]["alt_c"]
            graph.edges()[edge]["alt_c"] = None
            
            random.shuffle(free_ks)
            for ki, k in free_ks:
                O,D,_ = demands[k]
                try:
                    path1 = nx.dijkstra_path(graph,O,edge[0],"alt_c")
                    path2 = nx.dijkstra_path(graph,edge[1],D,"alt_c")
                except:
                    continue
                path = hf.nodes_to_edges_path(path1) + hf.nodes_to_edges_path(path2)
                epaths[ki] = path
                free_ks.remove((ki,k))
                break

            graph.edges()[edge]["alt_c"] = old_alt_c
        
        if len(with_edges) > 0:
            status = "not_solved"
            print("with_edges empty")
        else:
            print("with_edges not empty")
    
    
    X = epaths_to_X(graph,epaths,ks,demands)
 
    
    if status in ["neki","not_solved"]:
        zLD = (new_c@X.sum(axis=1))[0,0] - (vp["lam"].value @ vp["cap"])
    else: 
        zLD = -np.inf
    
    s = vp["cap"] - X.sum(axis=1).T
    
    ll = vp["lam"].value - s * float(alpha)
    ll = vp["lam"].value - s * float(alpha)
    ll = np.array(ll).flatten()
    ll[ll < 0] = 0 # lambda ne mora biti negativna
    if np.min(ll) < 0:
        pass
    vp["lam"].value = ll
    
    return (X, status, zLD, s)
