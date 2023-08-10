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



class Node:
    
    EPSI = 10e-3
    BETA = 1.1
    MAX_ITER_LR = 100
    MAX_UNDETERMINED_FEASIBILITIES = 5
    
    label = 0
    label_solved = 0
    
    
    def __init__(self, vp, graph, level = 0, parent = None, with_edges=[], without_edges=[], indetermined_feasibilities_left=0):
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
        
        self.with_edges = with_edges
        self.without_edges = without_edges
        
    
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
               "s":None,
               "message":""}
        
       
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
        
        
        
        # kopiram od starša, če:
            # -> sem with_edge otrok in smo vejali na uporabljenem edgu
            #( -> sem without_edge otrok in smo vejali po neuporabljenem edgu)
        
        if self.parent is not None and len(self.without_edges) == len(self.parent.without_edges):
            sol = dict(self.parent.sol)
            sol["status"] = "like_parent"
        else:
            q = 0 # initial number of iteration
            q_max = Node.MAX_ITER_LR # max number of iteration is q_max.
            while q <= q_max and beta > eps:
                ###############################
                # X, status, zLD, s = ss.dijkstra11(self.vp, graph, demands, alpha, self.with_edges, self.without_edges)
                # X, status, zLD, s = ss.dijkstra2(self.vp, graph, demands, self.U, self.L, alpha)
                # X, status, zLD, s = ss.cvxpy_linprog_LR(self.problem, self.vp, alpha)
                X, status, zLD, s = ss.dijkstra0_with_without(self.vp, graph, demands, alpha, self.with_edges, self.without_edges)
                
                ###############################################
                if status == "infeasible":
                    # print("conservation of flow constraint couldn't be satisfied at LD - infeasible")
                    self.sol = sol
                    self.sol["status"] = status
                    return
                
            
                values.append(zLD)
                
                if np.all(s >= 0) and status=="neki": # feasible for CLP
                    print("LR found feasible")
                    z = self.vp["c"].T @ np.sum(X,axis=1) # z
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
                            # cost too la
                            self.sol = sol
                            self.sol["status"] = "cost_to_large"
                            sol["zLD_ceil"] = max(int(np.ceil(LB)),int(LB_))
                            return
                        
                if flag > 2:
                    # print("beta halfed in ", q)
                    beta = beta/2
                    flag = 0
                
                # if zi == None:
                #     # print("zi is None")
                #     alpha = 1/(q+1)
                # else:
                #     alpha = abs(beta * (zi - zLD)/np.linalg.norm(s))
                alpha = 1/(q+1)
            
                q = q + 1
                
            print(q,beta)
            
            sol["zLD_ceil"] = max(np.ceil(LB),LB_)
            
            must_find_feasible = False
            
            if not must_find_feasible:
                sol["X"] = X if sol["X"] is None else sol["X"]
                sol["s"] = s
                self.sol = sol
                return values
            
            
            
            
            if status == "feasible" and sol["X"] is not None:
                sol["X"] = X if sol["X"] is None else sol["X"]
                sol["s"] = s
            
            elif self.indetermined_feasibilities_left > 0: 
                sol["status"] = "feasibility_unchecked"
            else:
                
                
                # poišči dopustno za CLP rešitev, s hevristiko ali eksaktno
                # (če ne obstaja je nedopustna)
                print("NODE WILL BE SEARCHED FOR FEASIBILITY")
                # sol["status"] += " NODE NOT SEARCHED FOR FEASIBILITY "
                
              
                self.problem_ex.solve()
                
                if self.problem_ex.status == "infeasible":
                    print("capacities cound't be satisfied - infeasible")
                    self.sol = sol
                    self.sol["status"] = self.problem_ex.status
                    return
                print(" exact success ")
                sol["status"] = "exact" # mislim, da tak nima otrok
                X = self.vp["X"].value
                
                
                z = int(self.vp["c"].T @ np.sum(X,axis=1))
                
                sol["cap_ok"] = True
                sol["X"] = X
                sol["z"] = z
            
        self.sol = sol
        
        return values
    
    
    
    


    
    def get_children(self):
        X = self.sol["X"]
        # s = self.sol["s"] if self.sol["s"] is not None else self.vp["cap"] - np.sum(X,axis=1)4
        
        # ce jih ze imamo
        if self.children is not None:
            return self.children
            
        
        used_edges_ids = list(np.where(np.sum(X,axis=1) > 0)[0])
        # naceloma tuki tko al tko ni without_edgov od parenta
        # etoi = lambda e: list(self.graph.edges()).index(e)
        itoe = lambda ei: list(self.graph.edges())[ei]
        candidate_edges = [itoe(ei) for ei in used_edges_ids if itoe(ei) not in self.with_edges]
        if len(candidate_edges) == 0:
            return None # TODO
        new_branching_edge = random.sample(candidate_edges,1)[0]
        
        
        ch1 = Node(self.vp, self.graph, level=self.level+1, parent=self,
                   with_edges=list(self.with_edges) + [new_branching_edge],
                   without_edges=list(self.without_edges),
                   indetermined_feasibilities_left=self.indetermined_feasibilities_left-1)
        ch2 = Node(self.vp, self.graph, level=self.level+1, parent=self,
                   with_edges=list(self.with_edges),
                   without_edges=list(self.without_edges) +[new_branching_edge],
                   indetermined_feasibilities_left=self.indetermined_feasibilities_left-1)

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
        stri += " left = " + str(self.indetermined_feasibilities_left)
        # stri += " " + repr(self.without_edges) + " "
        # if len(self.without_edges)>0:
        #     stri += repr(list(self.graph.edges())[self.without_edges[-1]])
        
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
    
    vp["lam"] = cp.Parameter(m, nonneg=True)
    

    # dolocim vrednosti parametrom
    vp["c"] = np.round(np.array([data["c"] for _,_, data in graph.edges(data=True)])) # BŠS
    vp["cap"] = np.floor(np.array([data["cap"] for _,_, data in graph.edges(data=True)])) # caps floor to int
    vp["B"] = -1 * np.array(nx.incidence_matrix(graph,oriented=True).todense())
    def demands_to_matrix(demands,n):
        H = np.zeros((n,len(demands)))
        for k,(Ok,Dk,d) in enumerate(demands):
            H[Ok,k] = d
            H[Dk,k] = -d
        return H      

    vp["H"] = demands_to_matrix(demands,n)
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
    UB = np.sum(vp["c"]) * vp["H"].shape[1] * np.max(vp["H"]) # to je vsi komoditiji grejo po vseh povezavah
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
        
        status = n.sol["status"]
        zLD_ceil = n.sol["zLD_ceil"]
        z = n.sol["z"]
        
        if status in ["infeasible","exact","cost_to_large"] or zLD_ceil >= UB:
            n.sol["message"] = "UB: {}".format(str(UB))
            continue
        
    
        
        if n.sol["cap_ok"]: # dopustna za prvotni CLP
            if z < UB:
                UB = z
                n_best = n
                
            n.sol["status"] = "feasible"
            if z == zLD_ceil:
                n.sol["status"] += " OPTIMAL for I"
                continue
        else:
            if n.sol["z"] == n.sol["zLD_ceil"]: # TODO sam pol se da popravit??
                n.sol["status"] += " OPTIMAL VALUE for I, NO SOLUTION"
                continue
            
        
        ch = n.get_children()
        if ch is not None:
            todo += ch
        else:
            break
        


    
        
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
