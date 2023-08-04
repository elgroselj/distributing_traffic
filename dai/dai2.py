import cvxpy as cp
from matplotlib import pyplot as plt
import numpy as np
import random
import networkx as nx
from scipy import sparse
from collections import Counter
from enum import Enum
from problem import Problem

import subproblem_solvers as ss

import importlib
importlib.reload(ss)



# class syntax

class Status(Enum):
    CLP_FEASIBLE = 1
    LR_FEASIBLE = 2
    INFEASIBLE = 3
    COST_TOO_LARGE = 4
    
class Sol:
    X = None
    zLD_ceil = None
    msg = ""
    status = None
    

class Node:
    
    EPSI = 10e-3
    BETA = 1.1
    MAX_ITER_LR = 100
    MAX_UNDETERMINED_FEASIBILITIES = 5
    
    label = 0
    label_solved = 0
    
    
    def __init__(self, vp, graph, level = 0, parent = None, without_edges=[], indetermined_feasibilities_left=0):
        self.vp = vp
        self.graph = graph
        
        self.sol = None
        self.children = None
        
        self.level = level
        self.label = Node.label
        Node.label += 1
        
        self.label_solved = None
        
        self.parent = parent
        
        self.indetermined_feasibilities_left = indetermined_feasibilities_left
        
        self.without_edges = without_edges
    
        _,m = self.vp["B"].value.shape
        _,t = self.vp["H"].value.shape
        
    
    def solve(self,LB_,UB_,graph,demands):
        # solve 
            # najde približek za zLD (najboljšo mejo, ki jo lahko da LA)
                # če ga ne najdemo s primarno hevristko
                # poskusimo eksaktno
                # možno je da je relaksacija nedopustna
            # vrne tudi najboljšo rešitev za CLP, na katero je naletel  
                # če je ni brne zadnjo, za pomoč pri vejanju
                
            
        self.label_solved = Node.label_solved
        Node.label_solved += 1
        
        # self.vp["lam"].value = np.zeros(self.vp["lam"].value.shape) # praktično nobene dopustne rešitve ne dobim med LR iteracijami
        
        sol = {"zLD_ceil":-np.inf,
               "status":"neki",
               "X":None,
               "cap_ok":False,
               "s":None}
        
       
        flag = 0 # povečamo ko smo v biližini LB => manjšati beto
        beta = Node.BETA #2 #TODO test
        eps = Node.EPSI
        alpha = 2
        
        # LB = 0  # če je zLD skos manjši od LB (ni izboljšanja) se beta zmanjšuje :(
        LB = -np.inf
        zi = None # vpliva na korak
        # zi = UB_ # meje so ful slabše tko -100
        sol["z"] = np.inf
        
        values = []
        q = 0 # initial number of iteration
        q_max = Node.MAX_ITER_LR # max number of iteration is q_max.
        
        # if self.branchingTF is not None and self.branchingTF[4] == False: # meja zDL mojega starša in dopustna rešitev mojega starša veljata zame
        if self.parent is not None and len(self.without_edges) == len(self.parent.without_edges):
            # print(" like_parent ")
            sol = dict(self.parent.sol)
            sol["status"] = " like_parent "
        else:
            while q <= q_max and beta > eps:
                ###############################
                X, status, zLD, s = ss.dijkstra11(self.vp, graph, demands, alpha, self.without_edges)
                # X, status, zLD, s = ss.dijkstra2(self.vp, graph, demands, self.U, self.L, alpha)
                # X, status, zLD, s = ss.cvxpy_linprog_LR(self.problem, self.vp, alpha)
                
                ###############################################
                if status == "infeasible":
                    # print("conservation of flow constraint couldn't be satisfied at LD - infeasible")
                    self.sol = sol
                    self.sol["status"] = status
                    return
                
            
                values.append(zLD)
                
                # if np.all(np.sum(X,axis=1) <= self.vp["cap"].value):#X* is feasible:
                if np.all(s >= 0):#X* is feasible for CLP an node:
                    print("LR found feasible")
                    z = self.vp["c"].value.T @ np.sum(X,axis=1) # z
                    zi = z
                    if z <= sol["z"]: # zapomnim si najcenejšo dopustno rešitev
                        sol["z"] = z
                        sol["X"] = X
                        sol["status"] = "feasible"
                        sol["cap_ok"] = True
                        sol["lam"] = self.vp["lam"].value
                        self.indetermined_feasibilities_left = Node.MAX_UNDETERMINED_FEASIBILITIES
                        
                        
                        

                if zLD < LB:
                    flag = 3
                else:
                    if zLD - LB < eps * max(1,LB):
                        flag = flag + 1
                    if zLD > LB:
                        LB = zLD
                        #beta = Node.BETA # TODO reset po izboljšanju
                        if LB > UB_:
                            # print(self.problem.status)
                            self.sol = sol
                            self.sol["status"] = "cost_to_large"
                            sol["zLD_ceil"] = max(int(np.ceil(LB)),int(LB_))
                            # print("cost_to_large")
                            # print("zLD:",zLD,"UB",UB_)
                            return
                        
                if flag > 2:
                    print("beta halfed in ", q)
                    beta = beta/2
                    flag = 0

                # s = self.vp["cap"].value - np.sum(X,axis=1)
                
                if zi == None:
                    # print("zi is None")
                    alpha /= 2
                else:
                    alpha = abs(beta * (zi - zLD)/np.linalg.norm(s))
                
                
                # print(alpha)
                
                # ll = self.vp["lam"].value - s * alpha
                # ll[ll < 0] = 0 # lambda ne mora biti negativna
                # self.vp["lam"].value = ll
                q = q + 1
                
            print(q,beta)
            
            sol["zLD_ceil"] = max(np.ceil(LB),LB_)
            
            must_find_feasible = False
            
            if not must_find_feasible:
                sol["X"] = X if sol["X"] is None else sol["X"]
                sol["s"] = s#self.vp["cap"].value - np.sum(sol["X"],axis=1)
                self.sol = sol
                return values
            
            
            
            
            if sol["X"] is not None:
                sol["X"] = X if sol["X"] is None else sol["X"]
                sol["s"] = s#self.vp["cap"].value - np.sum(sol["X"],axis=1)
                # sol["zLD_ceil"] = max(int(np.ceil(LB)),int(LB_))
                
                # vzamemo tesnejšo od mej (straši vs naša)
                
            elif self.parent is not None and len(self.without_edges) == len(self.parent.without_edges):#tj če nismo že prej prekinili
                zLD = sol["zLD"]
                sol = dict(self.parent.sol)
                sol["status"] += " like_parent "
                sol["zLD"] = zLD
            elif True or self.indetermined_feasibilities_left > 0: 
                sol["status"] = "feasibility_unchecked"
            else:
                
                
                # poišči dopustno za CLP rešitev, s hevristiko ali eksaktno
                # (če ne obstaja je nedopustna)
                print("NODE WILL BE SEARCHED FOR FEASIBILITY")
                # sol["status"] += " NODE NOT SEARCHED FOR FEASIBILITY "
                
                sub_problem = Problem
               
                
                
                def get_lower_upper_bounds_on_flow(k,lower):
                    l = {}
                    n = self
                    edges = list(graph.edges())
                    while n is not None and n.branchingTF is not None:
                        if n.branchingTF[1] == k and n.branchingTF[4] == lower:
                            ei = n.branchingTF[0]
                            e = edges[ei]
                            l[e] = n.branchingTF[3]
                        n = n.parent
                    return l
                
                
                def nodes_to_edges_path(inp):
                    nl = inp
                    el = []
                    for i in range(1,len(nl)):
                        el.append((nl[i-1], nl[i]))
                    return el
                
                # _,t = self.vp["H"].value.shape
                # _,m = self.vp["X"].value.shape
                
                
                
                
                
                X_dict = {}
                solution_found = True
                _,t = self.vp["H"].shape
                _,m = self.vp["B"].shape 
                
                remaining_cap = {e: graph.edges()[e]["cap"] for e in graph.edges()}
                for k in range(t):
                    Ok, Dk, num_k = demands[k]
                    # ubf
                    ubf = get_lower_upper_bounds_on_flow(k,lower=False)
                    remaining_cap_k = dict(remaining_cap)
                    for e in ubf:
                        remaining_cap_k[e] = min(remaining_cap[e], ubf[e])
                    nx.set_edge_attributes(graph, remaining_cap_k, "remaining_cap")
                
                    
                    # lbf
                    lbf = get_lower_upper_bounds_on_flow(k, lower=True)
                    
                    c = {e: graph.edges()[e]["c"] for e in graph.edges()}
                    for e in graph.edges():
                        if graph.edges()[e]["remaining_cap"] == 0:
                            c[e] = np.inf # nasičene / s kapaciteto 0 naj bodo neskončno drage
                    alt_c = {e: (0 if e in lbf else c[e]) for e in c} # te k rabijo bit uporabljene naj bojo zastonj
                    
                    nx.set_edge_attributes(graph, alt_c, "alt_c")
                    
                    On = list(graph.nodes())[Ok]
                    Dn = list(graph.nodes())[Dk]
                    
                    for ki in range(num_k):
                        
                        path_of_nodes = nx.dijkstra_path(graph,On,Dn,"alt_c")
                        path = nodes_to_edges_path(path_of_nodes)
                        length = sum([graph.edges()[e]["alt_c"] for e in path])
                        print(length)
                        etoei = lambda e : list(graph.edges()).index(e)
                        path_dict = {(etoei(e),k):1 for e in path}
                        X_dict = Counter(X_dict) + Counter(path_dict)
                        
                        for e in path:
                            if e in lbf:
                                lbf[e] -= 1
                                if lbf[e] == 0:
                                    del lbf[e]
                                    graph.edges()[e]["alt_c"] = alt_c[e]
                        
                        for e in path:
                            graph.edges()[e]["remaining_cap"] -= 1
                            remaining_cap[e] -= 1
                            if graph.edges()[e]["remaining_cap"] == 0:
                                graph.edges()[e]["alt_c"] = np.inf
                            elif remaining_cap[e] < 0:
                                # not solution
                                solution_found = False
                                break
                
                    if len(lbf) > 0:
                        # not solution    
                        print("lbf not empty")
                        solution_found = False
                        break
                             
                if solution_found:
                    print(" heuristic success")
                    sol["status"] += " heuristic "  
                    X = sparse.dok_matrix((m,t))
                    for a,k in X_dict:
                        X[a,k] = X_dict[(a,k)]
                    X = X.todense()
                else: # solve exactly
                    self.problem_ex.solve()
                    
                    if self.problem_ex.status == "infeasible":
                        print("capacities cound't be satisfied - infeasible")
                        self.sol = sol
                        self.sol["status"] = self.problem_ex.status
                        return
                    print(" exact success ")
                    sol["status"] += " exact "
                    X = self.vp["X"].value
                
                
                z = int(self.vp["c"].value.T @ np.sum(X,axis=1))
                
                sol["status"] += " feasible "
                sol["cap_ok"] = True
                sol["X"] = X
                sol["z"] = z
            
        self.sol = sol
        
        return values
    
    
    
    


    
    def get_children(self):
        X = self.sol["X"]
        # s = self.sol["s"] if self.sol["s"] is not None else self.vp["cap"].value - np.sum(X,axis=1)4
        
        # ce jih ze imamo
        if self.children is not None:
            return self.children
            
        
        used_edges = list(np.where(np.sum(X,axis=1) > 0)[0])
        # naceloma tuki tko al tko ni without_edgov od parenta
        new_without_edge = random.sample(used_edges,1)[0]
        
        ch1 = Node(self.vp, self.graph, level=self.level+1, parent=self, without_edges=list(self.without_edges),indetermined_feasibilities_left=self.indetermined_feasibilities_left-1)
        ch2 = Node(self.vp, self.graph, level=self.level+1, parent=self, without_edges=list(self.without_edges)+[new_without_edge],indetermined_feasibilities_left=self.indetermined_feasibilities_left-1)

        self.children = [ch1,ch2]
        
        return self.children
    
    def __str__(self, level=0):
        ret = "\t"*level+repr(self)+"\n"
        if self.children is not None:
            for child in self.children:
                ret += child.__str__(level+1)
        return ret
    
    def __repr__(self):
        # če še ni razvit
        if self.sol is None:
            # sprintaš samo kdaj je bil generiran
            return str(self.label) + ": not available"
        
        # sicer sprintaš tudi kdaj je bil razvit
        stri =  str(self.label) + " / " + str(self.label_solved) + ": " + self.sol["status"]
        
        # če smo ugotovili nedopustnost ja to to kar vemo
        if self.sol["status"] == "infeasible":
            return stri
        
        stri += " z: "+ str(self.sol["z"])+","
        stri += " zLD: " + str(self.sol["zLD_ceil"]) #+ "(" + str(self.sol["zLD_ceil_max"]) + ")"
        stri += " cap_ok = " + str(self.sol["cap_ok"])
        stri += " " + repr(self.without_edges) + " "
        if len(self.without_edges)>0:
            stri += repr(list(self.graph.edges())[self.without_edges[-1]])
        
        # če je razvejan, po čem
        # if self.children is not None:
            # a = self.children[0].branchingTF[0]
            # stri += " " + repr(list(self.graph.edges())[a]) + "(" + str(a) + ")"  + "," + repr(self.children[0].branchingTF[1:4]) + " "
            
            
        return stri

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
    
    vp["eta"] = cp.Parameter((m,t), nonneg = True)
    vp["xi"] = cp.Parameter((m,t), nonneg = True)

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
    
    return vp


def run(vp,graph,demands,MAX_ITER,MAX_ITER_LR):
    Node.label = 0
    Node.label_solved = 0
    Node.MAX_ITER_LR = MAX_ITER_LR
    n1 = Node(vp,graph)
    todo = [n1]

    n_best = None
    # celostevilska rešitev
    UB = np.sum(vp["c"].value) * vp["H"].value.shape[1] * np.max(vp["H"].value) # to je vsi komoditiji grejo po vseh povezavah
    q = 0
    while len(todo) > 0:
        if q >= MAX_ITER: break
        q += 1
        
        ############################################
            
        # if q % 2 == 0:
        #     index, _ = min(enumerate(todo), key=lambda tup : UB if tup[1].parent is None else tup[1].parent.sol["zLD_ceil"])
        # else:
        #     index = 0
        # n = todo.pop(index)
        
        #################################3
        # index, _ = max(enumerate(todo), key=lambda tup : len(tup[1].without_edges))
        # n = todo.pop(index)
        ######################################33
        n = todo.pop(0)
        ##############################################

        
        LB = n.parent.sol["zLD_ceil"] if n.parent is not None else 0
        # spodnjo mejo lahko vzameš od staršev, ker bo otrok imel kvečjemu dražje rešitve
        # LB DOBIŠ OD STARŠEV, UB DOBIŠ OD GLOBALNE CELE REŠITVE
        
        
        
        
        
        values = n.solve(LB,UB,graph,demands)
        
        
        try:
            plt.plot(values)
            # plt.show()
        except:
            pass
        
        
        if n.sol["status"] == "infeasible":
            continue
        
        
        zLD_ceil = n.sol["zLD_ceil"]
        z = n.sol["z"]
        
        if n.sol["status"] == "cost_to_large" or zLD_ceil >= UB: # ta veja bo samo še slabša (dražja) ali enaka
            n.sol["status"] += " COST too large("+str(UB)+")"
            continue
        
        if n.sol["cap_ok"]: # dopustna za prvotni CLP
            if z < UB:
                UB = z
                n_best = n
                
            n.sol["status"] += " FEASIBLE for I"
            if z == zLD_ceil:
                n.sol["status"] += " OPTIMAL for I"
                continue
        else:
            if n.sol["z"] == n.sol["zLD_ceil"]: # TODO sam pol se da popravit??
                n.sol["status"] += " OPTIMAL VALUE for I, NO SOLUTION"
                continue
            
        
        ch = n.get_children()
        todo += ch
        


    
        
    print(n1)
    
    if len(todo) == 0:
        print("VSE PREISKANO")
    else:
        print()
        print("Najmanjša cena ni manjša od: ",end="") # podati oceno koliko je še lufta do optimuma
        _, n = min(enumerate(todo), key=lambda tup : UB if tup[1].parent is None else tup[1].parent.sol["zLD_ceil"])
        end_LB = n.parent.sol["zLD_ceil"]
        print(end_LB)
        
    if n_best is not None:
        print(repr(n_best))
        print(sparse.csr_matrix(n_best.sol["X"]))
        

    return n_best, end_LB


def dai2_solve(graph,demands,MAX_ITER,MAX_ITER_LR):
    vp  = init_from_graph(graph,demands)

    n_best, end_LB = run(vp,graph,demands,MAX_ITER=MAX_ITER,MAX_ITER_LR=MAX_ITER_LR)
    
    message = "{} is LB".format(end_LB)
    if n_best is None:
        status = "no_feasible_solutions_found"
        X = None
        z = None
    else: 
        status, X, z = (n_best.sol["status"], n_best.sol["X"], n_best.sol["z"])
        
        
    return (status,X,z,message)
