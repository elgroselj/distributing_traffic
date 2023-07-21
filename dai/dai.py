# import cvxpy as cp
# from cvxopt import spmatrix
from matplotlib import pyplot as plt
import numpy as np
import random
import networkx as nx
from scipy import sparse
from collections import Counter
from docplex.mp.model import Model


    

class Node:
    
    EPSI = 10e-3
    BETA = 1.1
    MAX_ITER_LR = 100
    # TODO global lambda
    
    label = 0
    label_solved = 0
    
    def __init__(self, vp, graph, level = 0, parent = None, branchingTF = None):
        
        self.vp = vp
        self.graph = graph
        
        self.sol = None
        self.children = None
        
        self.level = level
        self.label = Node.label
        Node.label += 1
        
        self.label_solved = None
        
        self.parent = parent
        
        self.branchingTF = branchingTF
    
        

    
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
        
        # self.vp["lam"] = np.zeros(self.vp["lam"].shape) # praktično nobene dopustne rešitve ne dobim med LR iteracijami
        
        sol = {"zLD_ceil":-np.inf,
               "status":"neki",
               "X":None,
               "cap_ok":False,
               "s":None}
        
       
        flag = 0 # povečamo ko smo v biližini LB => manjšati beto
        beta = Node.BETA #2 #TODO test
        eps = Node.EPSI
        alpha = 1
        
        # LB = 0  # če je zLD skos manjši od LB (ni izboljšanja) se beta zmanjšuje :(
        LB = -np.inf
        zi = None # vpliva na korak
        # zi = UB_ # meje so ful slabše tko -100
        sol["z"] = np.inf
        
        values = []
        q = 0 # initial number of iteration
        q_max = Node.MAX_ITER_LR # max number of iteration is q_max.
        
        if self.branchingTF is not None and self.branchingTF[4] == False: # meja zDL mojega starša in dopustna rešitev mojega starša veljata zame
            # print(" like_parent ")
            sol = dict(self.parent.sol)
            sol["status"] = " like_parent "
        else:
            while q <= q_max and beta > eps:
                
                # self.problem.solve()
                X, zLD, status = cplex_solver(self.vp, self.branch_constraint_generator(), mode="LR")
                
                if status != "integer optimal solution":
                    print("put under infeasible: ", status)
        
                # if self.problem.status == "infeasible":
                    print("conservation of flow constraint couldn't be satisfied at LD - infeasible")
                    self.sol = sol
                    self.sol["status"] = "infeasible"
                    return
               
                
                    
                
                values.append(zLD)
                
                if np.all(np.sum(X,axis=1) <= self.vp["cap"]):#X* is feasible:

                    z = self.vp["c"].T @ np.sum(X,axis=1) # z
                    zi = z
                    if z <= sol["z"]: # zapomnim si najcenejšo dopustno rešitev
                        sol["z"] = z
                        sol["X"] = X
                        sol["status"] = "feasible"
                        sol["cap_ok"] = True
                        sol["lam"] = self.vp["lam"]
                        
                        # beta = 2 # TODO test # reset po izboljšanju
                        
                        
                        

                if zLD < LB:
                    flag = 3
                else:
                    if zLD - LB < eps * max(1,LB):
                        flag = flag + 1
                    if zLD > LB:
                        LB = zLD
                        # beta = Node.BETA # TODO reset po izboljšanju
                        # if LB > UB_:
                        #     # print(self.problem.status)
                        #     self.sol = sol
                        #     self.sol["status"] = "cost_to_large"
                        #     sol["zLD_ceil"] = max(int(np.ceil(LB)),int(LB_))
                        #     # print("cost_to_large")
                        #     # print("zLD:",zLD,"UB",UB_)
                        #     return
                        
                if flag > 2:
                    print("beta halfed in ", q)
                    beta = beta/2
                    flag = 0

                # s = self.vp["cap"] - np.sum(self.vp["X"],axis=1)
                s = self.vp["cap"] - np.sum(X,axis=1)
                
                if zi == None:
                    # print("zi is None")
                    alpha /= 2
                else:
                    alpha = abs(beta * (zi - zLD)/np.linalg.norm(s))
                
                
                # print(alpha)
                
                ll = self.vp["lam"] - s * alpha
                ll[ll < 0] = 0 # lambda ne mora biti negativna
                self.vp["lam"] = ll
                q = q + 1
                
            print(q,beta)
            
            # vzamemo tesnejšo od mej (straši vs naša)
            sol["zLD_ceil"] = max(np.ceil(LB),LB_)
            if sol["X"] is not None:
                sol["X"] = X if sol["X"] is None else sol["X"]
                sol["s"] = self.vp["cap"] - np.sum(sol["X"],axis=1)
                
                
            # elif self.branchingTF is not None and self.branchingTF[4] == False:#tj če nismo že prej prekinili
            #     zLD = sol["zLD"]
            #     sol = dict(self.parent.sol)
            #     sol["status"] += " like_parent "
            #     sol["zLD"] = zLD
            else:
                
                
                # TODO poišči dopustno za CLP rešitev, s hevristiko ali eksaktno
                # (če ne obstaja je nedopustna)
                # zaenkrat:
                print("NODE NOT SEARCHED FOR FEASIBILITY")
                # sol["status"] += " NODE NOT SEARCHED FOR FEASIBILITY "
                
                
               
                
                
                def get_lower_bounds_on_flow(k):
                    l = {}
                    n = self
                    edges = list(graph.edges())
                    while n is not None and n.branchingTF is not None:
                        if n.branchingTF[1] == k and n.branchingTF[4] == True:
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
                
                _,t = self.vp["H"].shape
                m = graph.number_of_edges()
                remaining_cap = {e: graph.edges()[e]["cap"] for e in graph.edges()}
                nx.set_edge_attributes(graph, remaining_cap, "remaining_cap")
                
                
                
                
                X_dict = {}
                solution_found = True
                for k in range(t):
                    Ok, Dk, num_k = demands[k]
                    lbf = get_lower_bounds_on_flow(k)
                    
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
                    # X = X.todense()
                else: # solve exactly
                    # self.problem_ex.solve()
                    X, _, status = cplex_solver(self.vp, self.branch_constraint_generator(), mode="feasibility")
                    
                    # if self.problem_ex.status == "infeasible":
                    if status != "integer optimal solution":
                        print("put under infeasible (exact-feasibility): ",status)
                        print("capacities cound't be satisfied - infeasible")
                        self.sol = sol
                        self.sol["status"] = "infeasible"
                        return
                    print(" exact success ")
                    sol["status"] += " exact "
                    # X = self.vp["X"]
                
                
                z = int(self.vp["c"].T @ np.sum(X,axis=1))
                
                sol["status"] += " feasible "
                sol["cap_ok"] = True
                sol["X"] = X
                sol["z"] = z
            
        self.sol = sol
        
        return values
    
    
    
    def branch_constraint_generator(self):
        n = self
        while n is not None and n.branchingTF is not None:
            yield n.branchingTF
            n = n.parent
        return True

    
    def get_children(self):
        X = self.sol["X"]
        s = self.sol["s"] if self.sol["s"] is not None else self.vp["cap"] - np.sum(X,axis=1)
        # ce jih ze imamo
        if self.children is not None:
            return self.children
        
        # preverimo, če so že vse celoštevilske
        # if np.all(np.abs(X - np.round(X)) < Node.EPSI): # vsi blizu celih
        #     self.children = []
        # else:
        
        # vejamo po spremenljivki, z največjim necelim delov PR 0.56 > 0,03
        # a,k = np.unravel_index(np.argmax(X % 1), )
        
        # random
        # a,k = np.unravel_index(np.random.randint(0,X.shape[0]*X.shape[1]),X.shape)
       
        # vejamo po spremenljivki glede na subgradient
        # a = np.argmax(s) # tam kjer je najmanj lufta sam ne uposteva dodanih omejitev TODO tkd, da na bit skos isto??
        # 1
        random.seed(123)
        # k = np.random.randint(0,X.shape[1])
        
        # 2 (hmm)
        # indeksi = list(np.where(X[a,:] > 0)[0])
        # k = random.sample(indeksi,1)[0] # omejimo eno od dobrin, ki uporablja to povezavo
        
        # 3
        # r = np.argsort(-s)
        # for a in r:
        #     k = np.random.randint(0,X.shape[1])
        #     val = int(X[a,k])
        #     branching = (a,k,val,val+1)
        #     if branching not in branches:
        #         self.branches.append(branching)
        #         break
        
        def is_new_branching(branching):
            n = self
            while n is not None and n.branchingTF is not None:
                if branching == n.branchingTF[:4]:
                    return False
                n = n.parent
            return True
        
        def create_branching(X,a,k):
            val = int(X[a,k])
            cap_a = int(self.vp["cap"][a])
            if cap_a <= val:
                val = cap_a - 1
            # tj če cap = 4, in mamo na njej val = 5: delimo <=3 >=4
            # tj če cap = 4, in mamo na njen val = 4: delimo <=3 >=4
            branching = (a,k,val,val+1)
            return branching
        
        branching = None
        if self.sol["cap_ok"]:
            # indeksi = list(np.where(X > 0)[0])
            indeksi = list(sparse.find(X > 0)[0])
            indeksi = random.sample(indeksi,len(indeksi))
            for i in indeksi:
                print("*",end="")
                a,k = np.unravel_index(i,X.shape)
                branching = create_branching(X,a,k)
                if is_new_branching(branching):
                    print()
                    break
        else:
            r = np.argsort(s)
            for a in r:
                k = np.random.randint(0,X.shape[1])
                branching = create_branching(X,a,k)
                if is_new_branching(branching):
                    break
        if branching is None:
            raise("branching is None")
        val = branching[2]
        
        # constraints1 = self.problem.constraints + [self.vp["X"][a,k] <= val] # false
        # constraints2 = self.problem.constraints + [self.vp["X"][a,k] >= (val +1)] # true        
        
        
        ch1 = Node(self.vp,self.graph,level=self.level+1,parent=self,branchingTF=branching + (False,))
        ch2 = Node(self.vp,self.graph,level=self.level+1,parent=self,branchingTF=branching + (True,))
        
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
        
        # če je razvejan, po čem
        if self.children is not None:
            a = self.children[0].branchingTF[0]
            stri += " " + repr(list(self.graph.edges())[a]) + "(" + str(a) + ")"  + "," + repr(self.children[0].branchingTF[1:4]) + " "
            
        return stri

def init_from_graph(graph,demands):
    # razberem dimenzije
    n = len(graph.nodes()) # 6 # |V|
    m = len(graph.edges()) # 10 # |E|
    t = len(demands) # 3
    
    vp = {} # slovar parametrov
    

    # dolocim vrednosti parametrom
    vp["c"] = np.round(np.array([data["c"] for _,_, data in graph.edges(data=True)])) # BŠS
    vp["cap"] = np.floor(np.array([data["cap"] for _,_, data in graph.edges(data=True)])) # caps floor to int
 
    vp["B"] = -1 * nx.incidence_matrix(graph,oriented=True)
    def demands_to_matrix_sparse(demands,n):
        H = sparse.dok_matrix((n,len(demands)))
        for k,(Ok,Dk,d) in enumerate(demands):
            H[Ok,k] = d
            H[Dk,k] = -d
        return H      
    vp["H"] = demands_to_matrix_sparse(demands,n)
    print(vp["H"])
    vp["lam"]= np.zeros(m)

    return vp





def run(vp,graph,demands,MAX_ITER,MAX_ITER_LR):
    Node.label = 0
    Node.label_solved = 0
    Node.MAX_ITER_LR = MAX_ITER_LR
    n1 = Node(vp,graph)
    L = [n1]

    n_best = None
    # celostevilska rešitev
    UB = np.sum(vp["c"]) * vp["H"].shape[1] * vp["H"].tocsr().max() # to je vsi komoditiji grejo po vseh povezavah
    q = 0
    while len(L) > 0:
        if q >= MAX_ITER: break
        q += 1
        
        ############################################
            
        if q % 2 == 0:
            index, _ = min(enumerate(L), key=lambda tup : UB if tup[1].parent is None else tup[1].parent.sol["zLD_ceil"])
        else:
            index = 0
        n = L.pop(index)
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
        L += ch
        


    
        
    print(n1)
    
    if len(L) == 0:
        print("VSE PREISKANO")
    else:
        print()
        print("Najmanjša cena ni manjša od: ",end="") # podati oceno koliko je še lufta do optimuma
        _, n = min(enumerate(L), key=lambda tup : UB if tup[1].parent is None else tup[1].parent.sol["zLD_ceil"])
        print(n.parent.sol["zLD_ceil"])
        
    if n_best is not None:
        print(repr(n_best))
        print(sparse.csr_matrix(n_best.sol["X"]))
        

    return n_best

  

def cplex_solver(vp, additional_constraints, mode):
    c = vp["c"]
    cap = vp["cap"]
    B = vp["B"]
    H = vp["H"]
    lam = vp["lam"]
    
    n,m = B.shape
    _,t = H.shape

    milp_model = Model(name = "myMILP")

    commodities = np.arange(t)

    edges = np.arange(m)#list(graph.edges())
    nodes = np.arange(n)#list(graph.nodes())

    X = milp_model.integer_var_matrix(edges, commodities, lb = 0, name="X")

    if mode == "exact":
        cost_expr = milp_model.sum(X[edge, com] * c[edge] for edge in edges for com in commodities)
    elif mode == "LR":
        cost_expr = (milp_model.sum(X[edge, com] * c[edge] for edge in edges for com in commodities) 
            #          +
            # milp_model.sum(lam[edge] * 
            #                (cap[edge] - 
            #                 milp_model.sum(
            #                     X[edge, com] for com in commodities)) for edge in edges)
            )
    elif mode == "feasibility":
        cost_expr = 0

    for node in nodes:
        for com in commodities:
            milp_model.add_constraint(milp_model.sum(X[edge, com] * B[node, edge] for edge in edges)  == H[node, com])


    if mode in ["exact","feasibility"]:
        for edge in edges:
            milp_model.add_constraint(milp_model.sum(X[edge, com] for com in commodities)  <=  cap[edge])
    

    for ad in additional_constraints:
        edge, com, val, valplus, side = ad
        if side == False:
            milp_model.add_constraint(X[edge, com] <= val)
        else:
            milp_model.add_constraint(X[edge, com] >= valplus)
            
        
    


    milp_model.set_objective("min", cost_expr)

    milp_model.print_information()
    
    
  
    
    sol = milp_model.solve()
    
    print(sol.solve_details)
    milp_model.print_solution()
    X_sparse = sparse.dok_matrix((m,t))
    for com in commodities:
        for edge in edges:
            if sol[X[edge,com]] != 0:
                X_sparse[edge,com] = sol[X[edge,com]]
    print(X_sparse)
    
    return (X_sparse, sol.objective_value, sol.solve_details.status) # TODO mogoče bo to da si ne zapominim vseh Xov fajn