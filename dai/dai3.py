import time
import cvxpy as cp
from matplotlib import cm, pyplot as plt
import numpy as np
import random
import networkx as nx
from scipy import sparse
from collections import Counter
from enum import Enum

import subproblem_solvers as ss
import helper_functions as hf

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
    
    EPSI = 10e-11
    BETA = 1.1
    MAX_ITER_LR = 100
    
    label = 0
    label_solved = 0
    
    
    # def __init__(self, obj, constraints, obj_ex, constraints_ex_additional, vp, graph, level = 0, parent = None, branchingTF = None):
    def __init__(self, vp, graph, level = 0, parent = None, branchingTF = None, U = None, L = None):
        self.vp = vp
        self.graph = graph
        
        self.sol = None
        self.children = None
        
        self.level = level
        self.label = Node.label
        Node.label += 1
        
        self.label_solved = None
        
        self.parent = parent
        
        
    
    def solve(self,LB_,UB_,graph,demands,timeout):
        # solve 
            # najde približek za zLD (najboljšo mejo, ki jo lahko da LA)
                # če ga ne najdemo s primarno hevristko
                # poskusimo eksaktno
                # možno je da je relaksacija nedopustna
            # vrne tudi najboljšo rešitev za CLP, na katero je naletel  
                # če je ni brne zadnjo, za pomoč pri vejanju
                
            
        self.label_solved = Node.label_solved
        Node.label_solved += 1
        
        # self.vp["lam"] = np.zeros(self.vp["lam"].shape) # praktično nobene dopustne rešitve ne dobim med LR iteracijami
        
        sol = {"zLD_ceil":-np.inf,
               "status":"neki",
               "X":None,
               "cap_ok":False,
               "s":None,
               "over_cap_count":None}
        
       
        flag = 0 # povečamo ko smo v biližini LB => manjšati beto
        beta = Node.BETA #2 #TODO test
        eps = Node.EPSI
        alpha = 2
        
        # LB = 0  # če je zLD skos manjši od LB (ni izboljšanja) se beta zmanjšuje :(
        LB = -np.inf
        zi = None # vpliva na korak
        # zi = UB_ # meje so ful slabše tko -100
        sol["z"] = np.inf
        min_over_cap_count = None
        min_over_cap_X = None
        min_over_cap_cost = None
        min_over_cap_node = None
        
        values = []
        lams = []
        alphas = []
        q = 0 # initial number of iteration
        q_max = Node.MAX_ITER_LR # max number of iteration is q_max.
        
        # iterations
        time_diff = 0
        before_t = time.perf_counter()
        LR_found_feasible = 0
        while q <= q_max and beta > eps and time_diff < timeout:
            time_diff = time.perf_counter() - before_t
            
            # print(time_diff)
            ###############################
            X, status, zLD, s = ss.dijkstra0(self.vp, graph, demands, alpha)
            # print(s[:,[0,6]])
            # print(s)
            # print(self.vp["lam"])           
            # hf.plot_solution_graph(self.graph,X,edge_label1="alt_c",font_size=12,demands=demands,node_size=150,figure_size=(5,5))
            # X, status, zLD, s = ss.astar(self.vp, graph, demands, alpha)
            lams.append(np.copy(self.vp["lam"]))
            ###############################################
            if status == "infeasible":
                # print("conservation of flow constraint couldn't be satisfied at LD - infeasible")
                self.sol = sol
                self.sol["status"] = status
                return
            
            
            values.append(zLD)
            
            # if np.all(np.sum(X,axis=1) <= self.vp["cap"]):#X* is feasible:
            if np.all(s >= 0):#X* is feasible for CLP an node:
                LR_found_feasible += 1
                print("LR found feasible")
                z = float(self.vp["c"].T @ X.sum(axis=1))
                zi = z
                if z <= sol["z"]: # zapomnim si najcenejšo dopustno rešitev
                    sol["z"] = z
                    sol["X"] = X
                    sol["status"] = "feasible"
                    sol["cap_ok"] = True
                    sol["lam"] = self.vp["lam"]
                    sol["over_cap_count"] = 0
            elif zi is None:
                over_cap_count = -np.sum(s[s<0])
                over_cap_cost = int(self.vp["c"].T @ X.sum(axis=1))
                if (min_over_cap_count is None or over_cap_count < min_over_cap_count or 
                    (over_cap_count == min_over_cap_count and over_cap_cost < min_over_cap_cost)):
                    min_over_cap_count = over_cap_count
                    min_over_cap_cost = over_cap_cost
                    min_over_cap_X = X
                
            if zLD < LB:
                flag = 3
            else:
                if zLD - LB < eps * max(1,LB):
                    flag = flag + 1
                if zLD > LB:
                    LB = zLD
                    
                    
            if flag > 2:
                print("beta halfed in ", q)
                if zi is not None:
                    beta = beta/2
                flag = 0

            
            if zi == None:
                # alpha /= 2
                alpha = 1/(q+1)
            else:
                alpha = abs(beta * (zi - zLD)/np.linalg.norm(s))
            alphas.append(alpha)
            
            
            q = q + 1
                
        print(q,beta)
        print("LR_found_feasible",LR_found_feasible)
        
        sol["zLD_ceil"] = max(np.ceil(LB),LB_)
        if sol["X"] is None:
            
            
            # TODO poišči dopustno za CLP rešitev, s hevristiko ali eksaktno
            # (če ne obstaja je nedopustna)
            # zaenkrat:
            print("NODE NOT SEARCHED FOR FEASIBILITY")
            # sol["status"] += " NODE NOT SEARCHED FOR FEASIBILITY "
            
            
            ks = []
            for ki, (_,_,num_k) in enumerate(demands):
                for i in range(num_k):
                    ks.append(ki)
            X = hf.assemble(ks,graph,demands)
            
            if X is not None:
                z = int(self.vp["c"].T @ np.sum(X,axis=1))
                sol["z"] = z
            
            # sol["status"] += " feasible "
                sol["status"] = "feasible"
                
                sol["cap_ok"] = True
                sol["X"] = X
            else:
                if min_over_cap_X is None:
                    sol["status"] = "no_feasible_solution_found"
                else:
                    sol["status"] = "over_cap"
                    sol["X"] = min_over_cap_X
                    sol["z"] = min_over_cap_cost
                    sol["over_cap_count"] = min_over_cap_count
            
            
        self.sol = sol
        
        return values,lams,alphas
    

    
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
        
        # če je razvejan, po čem
        if self.children is not None:
            a = self.children[0].branchingTF[0]
            stri += " " + repr(list(self.graph.edges())[a]) + "(" + str(a) + ")"  + "," + repr(self.children[0].branchingTF[1:4]) + " "
            
        return stri

def init_from_graph(graph,demands,lam0=None):
    # razberem dimenzije
    n = len(graph.nodes()) # 6 # |V|
    m = len(graph.edges()) # 10 # |E|
    t = len(demands) # 3
    
    vp = {
        "X":None,
        "c":None,
        "cap":None,
        "B":None,
        "H":None,
        "lam":None,
          } # slovar spremenljivk in parametrov


    # dolocim vrednosti parametrom
    vp["c"] = np.round(np.array([data["c"] for _,_, data in graph.edges(data=True)])) # BŠS
    vp["cap"] = np.floor(np.array([data["cap"] for _,_, data in graph.edges(data=True)])) # caps floor to int
    vp["B"] = -1 * nx.incidence_matrix(graph,oriented=True)
    def demands_to_matrix(demands,n):
        H = sparse.dok_matrix((n,len(demands)))
        for k,(Ok,Dk,d) in enumerate(demands):
            H[Ok,k] = d
            H[Dk,k] = -d
        return H.tolil()      
    vp["H"] = demands_to_matrix(demands,n)
    vp["lam"] = sparse.lil_array(np.zeros((m))) if lam0 is None else sparse.lil_array(lam0)


    return vp


def run(vp,graph,demands,MAX_ITER,MAX_ITER_LR,timeout):
    Node.label = 0
    Node.label_solved = 0
    Node.MAX_ITER_LR = MAX_ITER_LR
    n = Node(vp,graph)

    n_best = None
    # celostevilska rešitev
    UB = np.sum(vp["c"]) * vp["H"].shape[1] * vp["H"].tocsr().max()# to je vsi komoditiji grejo po vseh povezavah
    
    
    
    

    
    LB = 0
    # spodnjo mejo lahko vzameš od staršev, ker bo otrok imel kvečjemu dražje rešitve
    # LB DOBIŠ OD STARŠEV, UB DOBIŠ OD GLOBALNE CELE REŠITVE
    
    
    
    
    UB = np.sum(vp["c"]) * vp["H"].tocsr().max(axis=0).sum() # to je vsi komoditiji grejo po vseh povezavah
    returned = n.solve(LB,UB,graph,demands,timeout)
    
    print("optimal lambda")
    print(vp["lam"])
    
    if False:
        plt.rcParams.update({'font.size': 17})
        try:
            values, lams, alphas = returned
        except:
            pass
        
        try:
            plt.title("konvergenca k v*")
            plt.xlabel("iteracija")
            plt.ylabel("približek za v*")
            plt.plot(values)
            plt.show()
        except:
            pass
        
        try:
            plt.title("konvergenca lambde po komponentah")
            plt.xlabel("iteracija")
            plt.ylabel("vrednosti komponent lambde")
            plt.plot([lam[0] for lam in lams])
            plt.show()
        except:
            pass
        
        # try:
        #     plt.plot(alphas)
        #     plt.show()
        # except:
        #     pass
        
        # try:
        #     for i in range(len(lams[0][0,:])):
        #         plt.scatter([lam[:,i] for lam in lams],values,c=["k"]*len(values))
        #         plt.show()
        # except:
        #     pass
        
        
        norm_values = np.array(values)
        norm_values = norm_values- min(norm_values)
        norm_values = norm_values/ max(norm_values)
        x = [lam[:,0] for lam in lams]
        y = [lam[:,6] for lam in lams]
        plt.scatter(x,y,c=values)
        # plt.scatter(x,y,s=norm_values*20)
        for i,(xi, yi) in enumerate(zip(x, y)):
            if i > 13 : i = None
            plt.text(xi, yi, i, va='bottom', ha='center')
        
        plt.xlabel("lambda_1")
        plt.ylabel("lambda_7")
        plt.title("konvergenca k v*")
        cbar = plt.colorbar()
        cbar.set_label('približek za v*')
        plt.show()
        
        
        
    
    
    
    zLD_ceil = n.sol["zLD_ceil"]
    z = n.sol["z"]
    
    X = n.sol["X"]
    # s = vp["cap"] - X.sum(axis=1).T
    # if not np.all(s >= 0):
    #     print("false_feasible")
    
    
    if n.sol["cap_ok"]: # dopustna za prvotni CLP
        if z < UB:
            UB = z
            n_best = n
            
        # n.sol["status"] += " FEASIBLE for I"
        n.sol["status"] = "feasible"
        if z == zLD_ceil:
            # n.sol["status"] += " OPTIMAL for I"
            n.sol["status"] = "optimal"
    # else:
        # if n.sol["z"] == n.sol["zLD_ceil"]: # TODO sam pol se da popravit??
        #     # n.sol["status"] += " OPTIMAL VALUE for I, NO SOLUTION"
        #     n.sol["status"] = "optimal_but_no_solution_found"
        # if n.sol["status"] == "over_cap":
        #     n.sol["X"]
        


    
        

    print("Najmanjša cena ni manjša od: ",end="") # podati oceno koliko je še lufta do optimuma
    end_LB = n.sol["zLD_ceil"]
    print(end_LB)
        
    if n_best is not None:
        print(repr(n_best))
        print(sparse.csr_matrix(n_best.sol["X"]))   

    return n, end_LB, 


def dai3_solve(graph,demands,MAX_ITER,MAX_ITER_LR,timeout,lam0 = None):
    vp  = init_from_graph(graph,demands,lam0=lam0)

    n_best, end_LB = run(vp,graph,demands,MAX_ITER=MAX_ITER,MAX_ITER_LR=MAX_ITER_LR,timeout=timeout)
    
    message = ""
    LB = end_LB
    if n_best is None:
        status = "no_feas"
        X = None
        z = None
        over_cap_count = None
    else: 
        status, X, z = (n_best.sol["status"], n_best.sol["X"], n_best.sol["z"])
        over_cap_count = n_best.sol["over_cap_count"]
        
        
    return (status,X,z,message,LB,over_cap_count)
