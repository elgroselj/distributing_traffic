import cvxpy as cp
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from node import Node 
COLORS = "brgpy"

def plot_multigraph(graph, with_labels=True,font_size=5,figure_size=(20,20)):
    plt.figure(figsize=figure_size)
    G = nx.MultiDiGraph(graph)
    try:
        pos = {n:(G.nodes[n]["x"],G.nodes[n]["y"]) for n in G.nodes()}
    except:
        pos = nx.circular_layout(G)
    nx.draw_networkx_nodes(G, pos, node_size=0, alpha=1)
    nx.draw_networkx_labels(G, pos, font_size=font_size)
    ax = plt.gca()
    for e in G.edges:
        edge_color = G.edges[e]["color"] if "color" in G.edges[e] else "k"
        ax.annotate("",
                    xy=pos[e[1]], xycoords='data',
                    xytext=pos[e[0]], textcoords='data',
                    arrowprops=dict(arrowstyle="->", color=edge_color,
                                    shrinkA=5, shrinkB=5,
                                    patchA=None, patchB=None,
                                    connectionstyle="arc3,rad=rrr".replace('rrr', str(0.1 * (e[2] + 1))
                                                                           ),
                                    ),
                    )
    if with_labels:
        digraph = nx.DiGraph(G)
        nx.draw_networkx_edge_labels(digraph, pos,
                                     edge_labels={e: "{},{}".format(int(np.round(digraph.edges[e]["c"])), int(np.floor(digraph.edges[e]["cap"])))
                                                  for e in list(digraph.edges)},
                                     font_size=font_size)

    plt.axis('off')
    plt.show()

def plot_solution_graph(graph,X,with_labels=True,font_size=5,figure_size=(20,20)):
    multi = nx.MultiDiGraph(graph)
    print("k\tCOLOR")
    for k in range(X.shape[1]):
        path = [e for i,e in enumerate(graph.edges) if X[i,k] != 0]
        multi.add_edges_from(path, color=COLORS[k])
        print(k,"\t",COLORS[k])
    # print(list(multi.edges()))
    plot_multigraph(multi,with_labels=with_labels,font_size=font_size,figure_size=figure_size)

def init_from_graph(graph,demands):
    # razberem dimenzije
    n = len(graph.nodes()) # 6 # |V|
    m = len(graph.edges()) # 10 # |E|
    t = len(demands) # 3
    
    vp = {} # slovar spremenljivk in parametrov
    
    # spremenljivke
    vp["X"] = cp.Variable((m,t),integer=True)

    # parametri
    vp["c"] = cp.Parameter(m, integer=True)
    vp["cap"] = cp.Parameter(m, integer = True)
    vp["B"] = cp.Parameter((n,m), integer = True)
    vp["H"] = cp.Parameter((n,t), integer= True)
    vp["lam"] = cp.Parameter(m, nonneg=True)

    # dolocim vrednosti parametrom
    vp["c"].value = np.round(np.array([data["c"] for _,_, data in graph.edges(data=True)])) # BŠS
    vp["cap"].value = np.floor(np.array([data["cap"] for _,_, data in graph.edges(data=True)])) # caps floor to int
    print(vp["c"].value)
    print(vp["cap"].value)
    vp["B"].value = -1 * np.array(nx.incidence_matrix(graph,oriented=True).todense())
    def demands_to_matrix(demands,n):
        H = np.zeros((n,len(demands)))
        for k,(Ok,Dk,d) in enumerate(demands):
            H[Ok,k] = d
            H[Dk,k] = -d
        return H      
    # vp["H"].value = np.array([[1, 3, 2],
    #                     [0, 0, 0],
    #                     [0, 0, 0],
    #                     [-1, 0, 0],
    #                     [0, -3, 0],
    #                     [0, 0, -2]
    #                     ])
    vp["H"].value = demands_to_matrix(demands,n)
    print(vp["H"].value)
    vp["lam"].value = np.zeros(m)



    # kriterijska funkcija
    obj = cp.Minimize(vp["c"].T @ cp.sum(vp["X"],axis=1) + vp["lam"] @ (cp.sum(vp["X"],axis=1) - vp["cap"]))

    # omejitve
    constraints = [
        vp["B"] @ vp["X"] == vp["H"],
        vp["X"] >= 0   
    ]

    # prob = cp.Problem(obj, constraints)
    # # print(prob.is_dpp()) TODO
    return (obj, constraints, vp)


def run(obj,constraints,vp,graph,MAX_ITER,INIT_NUM_STEPS):
    Node.label = 0
    Node.label_solved = 0
    Node.INIT_NUM_STEPS = INIT_NUM_STEPS
    n1 = Node(obj,constraints,vp)
    L = [n1]

    n_best = None
    UB = np.inf # celostevilski
    
    q = 0
    while len(L) > 0:
        if q >= MAX_ITER: break
        q += 1
        
        # 1
        # prednost imajo otroci staršev z cap_ok == True
        n = None
        for n2 in L:
            if n2.parent is not None and "cap_ok" in n2.parent.sol and n2.parent.sol["cap_ok"]:
                n = n2
                L.remove(n2)
                break
        if n is None:
            n = L.pop(0)
        
        # 2
        # sortiramo najprej po cap_ok, potem pa kdo ima najbolj ohlapno mejo
        # def f(n_):
        #     if n_.parent is None:
        #         return (0,UB) 
        #     else:
        #         return (~n_.parent.sol["cap_ok"], n_.parent.sol["zLD"]) # prvo (0 == True, majhna št)
        # L.sort(key=lambda n_ : f(n_)) 
        # n = L.pop(0)

        values = n.solve(UB)
        try:
            plt.plot(values)
        except:
            pass
        
        # print(n.sol["status"])
        if n.sol["status"] == "infeasible":
            continue
        
        # zLD = n.sol["zLD"] # optimistična hevristika
        zLD_ceil = n.sol["zLD_ceil"] # zaokrožimo BŠS, dobimo bolj tesno mejo
        
        z = n.sol["z"]
        
        if z > UB: # ta veja bo samo še slabša (dražja)
            n.sol["status"] += " COST too large("+str(UB)+")"
            continue
        
        if n.sol["cap_ok"]: # dopustna za prvotni CLP
            if z < UB:
                UB = z
                # X_best = n.sol["X"]
                n_best = n
                # LB_best = zLD_ceil
            n.sol["status"] += " FEASIBLE for I"
            if z == zLD_ceil:
                n.sol["status"] += " OPTIMAL for I"
                continue
        
        ch = n.get_children()
        L += ch
        


    if len(L) == 0: print("VSE PREISKANO")
    print(n1)
    if n_best is not None:
        print(repr(n_best))
        print(repr(n_best.sol["X"]))
        
    # try:        
        
    #     print(n_best.sol["status"],": ", UB)
    #     print(repr(X_best))
    #     print("LB_best: ", LB_best)
    #     print(np.sum(X_best,axis=1) <= vp["cap"].value)
    #     print(vp["c"].value.T @ np.sum(X_best,axis=1))
    # except: pass


    # TODO bug ko "vse preiskano" pa ni optimalne rešitve
    return n_best # TODO podati oceno koliko je še lufta do optimuma

def run2(obj,constraints,vp,MAX_ITER):
    q = 0 # initial number of iteration
    flag = 0
    betha = 2
    q_max = MAX_ITER # max number of iteration is q_max.
    UB = np.sum(vp["c"].value) * vp["H"].value.shape[1] * np.max(vp["H"].value) # to je vsi komoditiji grejo po vseh povezavah
    LB = -np.inf  # initial upper bound and lower bound 
    eps = 10**(-11)
    
    UB_min = UB
    X_best = None
    
    problem = cp.Problem(obj, constraints)
    
    # Node.label = 0
    # Node.INIT_NUM_STEPS = INIT_NUM_STEPS
    # n = Node(obj,constraints,vp)
    values = []
    while q <= q_max and betha > eps:
        problem.solve()
        X = vp["X"].value
        zLD = problem.value
        values.append(zLD)
        if np.all(np.sum(X,axis=1) <= vp["cap"].value):#X* is feasible:
            UB = vp["c"].value.T @ np.sum(X,axis=1) # z
            if UB < UB_min:
                UB_min = UB
                X_best = X
                
        if zLD < LB:
            flag = 3
        else:
            if zLD - LB < eps * max(1,LB):
                flag = flag + 1
            if zLD > LB:
                LB = zLD
                
        if flag > 2:
            betha = betha/2
            flag = 0

        s = vp["cap"].value - np.sum(vp["X"].value,axis=1)
        alpha = betha * (UB - zLD)/np.linalg.norm(s)
        ll = vp["lam"].value - s * alpha
        ll[ll < 0] = 0 # lambda ne mora biti negativna
        vp["lam"].value = ll
        q = q + 1
    print(q,betha)
    try:
        plt.plot(values)
    except:
        pass
    return (LB,UB,X_best)
    

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
