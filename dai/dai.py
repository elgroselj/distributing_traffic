import cvxpy as cp
from matplotlib import pyplot as plt
import numpy as np
import random
import networkx as nx
from scipy import sparse

class Node:
    
    EPSI = 10e-11
    BETA = 1.1
    MAX_ITER_LR = 100
    
    label = 0
    label_solved = 0
    
    def __init__(self, obj, constraints, vp, graph, level = 0, parent = None, branchingTF = None):
        self.problem = cp.Problem(obj, constraints)
        self.vp = vp
        self.graph = graph
        
        self.sol = None
        self.children = None
        
        self.level = level
        self.label = Node.label
        Node.label += 1
        
        self.label_solved = None
        
        self.parent = parent
        
        # self.branches = branches
        self.branchingTF = branchingTF
    
        

    
    def solve(self,LB_,UB_):
        # solve 
            # najde približek za zLD (najboljšo mejo, ki jo lahko da LA)
                # če ga ne najdemo s primarno hevristko
                # poskusimo eksaktno
                # možno je da je relaksacija nedopustna
            # vrne tudi najboljšo rešitev za CLP, na katero je naletel  
                # če je ni brne zadnjo, za pomoč pri vejanju
            
        self.label_solved = Node.label_solved
        Node.label_solved += 1
        
        # self.vp["lam"].value = np.zeros(self.vp["lam"].value.shape)
        
        sol = {"zLD_ceil":-np.inf,
               "status":"neki",
               "X":None,
               "cap_ok":False,
               "solved_exacly":False}
        
       
        flag = 0 # povečamo ko smo v biližini LB => manjšati beto
        beta = Node.BETA #2 #TODO test
        eps = Node.EPSI
        
        LB = 0
        zi = None # vliva na korak
        sol["z"] = np.inf
        
        values = []
        q = 0 # initial number of iteration
        q_max = Node.MAX_ITER_LR # max number of iteration is q_max.
        
        while q <= q_max and beta > eps:
            
            self.problem.solve()
            
    
            if self.problem.status == "infeasible":
                print("conservation of flow constraint couldn't be satisfied at LD - infeasible")
                self.sol = sol
                self.sol["status"] = self.problem.status
                return
            
            X = self.vp["X"].value
            zLD = self.problem.value
            
            if zLD > UB_:
                print(self.problem.status)
                self.sol = sol
                self.sol["status"] = "cost_to_large"
                print("cost_to_large")
                print("zLD:",zLD,"UB",UB_)
                break
            
            values.append(zLD)
            
            if np.all(np.sum(X,axis=1) <= self.vp["cap"].value):#X* is feasible:

                z = self.vp["c"].value.T @ np.sum(X,axis=1) # z
                zi = z
                if z <= sol["z"]:
                    sol["z"] = z
                    sol["X"] = X
                    sol["status"] = "feasible"
                    sol["cap_ok"] = True
                    
                    # beta = 2 # TODO test # reset po izboljšanju
                    
                    
                    
                    
            if zLD < LB:
                flag = 3
            else:
                if zLD - LB < eps * max(1,LB):
                    flag = flag + 1
                if zLD > LB:
                    LB = zLD
                    
            if flag > 2:
                # print("beta halfed")
                beta = beta/2
                flag = 0

            s = self.vp["cap"].value - np.sum(self.vp["X"].value,axis=1)
            
            if zi == None:
                print("zi is None")
                alpha = 0.005
            else:
                alpha = abs(beta * (zi - zLD)/np.linalg.norm(s))
            
            
            # print(alpha)
            
            ll = self.vp["lam"].value - s * alpha
            ll[ll < 0] = 0 # lambda ne mora biti negativna
            self.vp["lam"].value = ll
            q = q + 1
            
        print(q,beta)
        
        sol["X"] = X if sol["X"] is None else sol["X"]
        sol["s"] = self.vp["cap"].value - np.sum(sol["X"],axis=1)
        sol["zLD_ceil"] = max(int(np.ceil(LB)),int(LB_))
        # vzamemo tesnejšo od mej (straši vs naša)
        
        self.sol = sol
        
        return values
    
    
    
    


    
    def get_children(self):
        X = self.sol["X"]
        s = self.sol["s"]
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
            cap_a = int(self.vp["cap"].value[a])
            if cap_a <= val:
                val = cap_a - 1
            # tj če cap = 4, in mamo na njej val = 5: delimo <=3 >=4
            # tj če cap = 4, in mamo na njen val = 4: delimo <=3 >=4
            branching = (a,k,val,val+1)
            return branching
        
        branching = None
        if self.sol["cap_ok"]:
            indeksi = list(np.where(X > 0)[0])
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
        
        constraints1 = self.problem.constraints + [self.vp["X"][a,k] <= val] # false
        constraints2 = self.problem.constraints + [self.vp["X"][a,k] >= (val +1)] # true        
        
        
        ch1 = Node(self.problem.objective,constraints1,self.vp,self.graph,level=self.level+1,parent=self,branchingTF=branching + (False,))
        ch2 = Node(self.problem.objective,constraints2,self.vp,self.graph,level=self.level+1,parent=self,branchingTF=branching + (True,))
        
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


def run(obj,constraints,vp,graph,MAX_ITER,MAX_ITER_LR):
    Node.label = 0
    Node.label_solved = 0
    n1 = Node(obj,constraints,vp,graph)
    L = [n1]

    n_best = None
    # celostevilska rešitev
    UB = np.sum(vp["c"].value) * vp["H"].value.shape[1] * np.max(vp["H"].value) # to je vsi komoditiji grejo po vseh povezavah
    q = 0
    while len(L) > 0:
        if q >= MAX_ITER: break
        q += 1
        
        ############################################
        def f(n_):
            if n_.parent is None:
                return (not True,UB) 
            else:
                print("*",end="")
                # print(~n_.parent.sol["cap_ok"])
                # return (~n_.parent.sol["cap_ok"], n_.parent.sol["zLD_ceil"]) # prvo (0 == True, majhna št)
                return (not n_.parent.sol["cap_ok"], n_.parent.sol["zLD_ceil"]) # prvo (0 == True, majhna št)
            
        if q % 2 == 0:
            L.sort(key=lambda n_ : f(n_))
        else:
            L.sort(key=lambda n_ : n_.label)
        print([(f(n_), n_.label) for n_ in L])
        n = L.pop(0)
        ##############################################

        
        LB = n.parent.sol["zLD_ceil"] if n.parent is not None else 0
        # spodnjo mejo lahko vzameš od staršev, ker bo otrok imel kvečjemu dražje rešitve
        # LB DOBIŠ OD STARŠEV, UB DOBIŠ OD GLOBALNE CELE REŠITVE
        
        
        Node.MAX_ITER_LR = MAX_ITER_LR
        
        
        values = n.solve(LB,UB)
        
        
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
        L.sort(key=lambda n_ : f(n_))
        print()
        print("Najmanjša cena ni manjša od: ",end="") # podati oceno koliko je še lufta do optimuma
        print(L.pop(0).parent.sol["zLD_ceil"])
        
    if n_best is not None:
        print(repr(n_best))
        print(sparse.csr_matrix(n_best.sol["X"]))
        

    return n_best

# def run2(obj,constraints,vp,MAX_ITER):
#     q = 0 # initial number of iteration
#     flag = 0
#     beta = 2
#     q_max = MAX_ITER # max number of iteration is q_max.
#     UB = np.sum(vp["c"].value) * vp["H"].value.shape[1] * np.max(vp["H"].value) # to je vsi komoditiji grejo po vseh povezavah
#     LB = -np.inf  # initial upper bound and lower bound 
#     eps = 10**(-11)
    
#     UB_min = UB
#     X_best = None
    
#     problem = cp.Problem(obj, constraints)
    
#     # Node.label = 0
#     # Node.INIT_NUM_STEPS = INIT_NUM_STEPS
#     # n = Node(obj,constraints,vp)
#     values = []
#     while q <= q_max and beta > eps:
#         problem.solve()
#         X = vp["X"].value
#         zLD = problem.value
#         values.append(zLD)
#         if np.all(np.sum(X,axis=1) <= vp["cap"].value):#X* is feasible:
#             UB = vp["c"].value.T @ np.sum(X,axis=1) # z
#             if UB < UB_min:
#                 UB_min = UB
#                 X_best = X
                
#         if zLD < LB:
#             flag = 3
#         else:
#             if zLD - LB < eps * max(1,LB):
#                 flag = flag + 1
#             if zLD > LB:
#                 LB = zLD
                
#         if flag > 2:
#             beta = beta/2
#             flag = 0

#         s = vp["cap"].value - np.sum(vp["X"].value,axis=1)
#         alpha = beta * (UB - zLD)/np.linalg.norm(s)
#         ll = vp["lam"].value - s * alpha
#         ll[ll < 0] = 0 # lambda ne mora biti negativna
#         vp["lam"].value = ll
#         q = q + 1
#     print(q,beta)
#     try:
#         plt.plot(values)
#     except:
#         pass
#     return (LB,UB,X_best)