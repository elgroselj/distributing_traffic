import cvxpy as cp
from matplotlib import pyplot as plt
import numpy as np
import random

class Node:
    EPSI = 10e-11
    INIT_NUM_STEPS = 100
    label = 0
    label_solved = 0
    def __init__(self, obj, constraints, vp, level = 0, parent = None, branchingTF = None):
        self.problem = cp.Problem(obj, constraints)
        self.vp = vp
        self.sol = None
        self.children = None
        
        self.level = level
        self.label = Node.label
        Node.label += 1
        
        self.label_solved = None
        
        self.parent = parent
        
        # self.branches = branches
        self.branchingTF = branchingTF
    
        
    
    def solve(self,UB,lam0 = None):
        self.label_solved = Node.label_solved
        Node.label_solved += 1
        
        sol = {"zLD":-np.inf,"status":None,"X":None,"lam":None,"cap_ok":False}
        no_change_counter = 0
        betha = 2
        if lam0 is None:
            self.vp["lam"].value = np.zeros(self.vp["lam"].value.shape)
        else:
            self.vp["lam"].value = lam0
        
        tt_init = Node.INIT_NUM_STEPS #100 # pomenbno je, da dobro skonvergiramo na začetku TODO analiziraj obnašanje
        tt = max(20,int(tt_init/(self.level+1)))
        print(tt)
        
        values = []
        for t in range(tt):
            # print(t,":")
            # rešimo LD pri neki lambdi
            # TODO reši LD bolj učinkovito     
            self.problem.solve()
            
            # če LD slučajno ni rešljiv, bomo šli ven že v 1. koraku
            # to je ko pogojem o ohranitvi toka ne more bit zadoščeno
            if self.problem.status == "infeasible":
                print("conservation of flow constraint couldn't be satisfied at LD - infeasible")
                self.sol = sol
                self.sol["status"] = self.problem.status
                return
            
            if self.problem.status == "unbounded":
                raise("LD - unbounded")
                break
            
            # zapomnimo si najtesnejšo (najvišjo) rešitev LD
            if self.problem.value > sol["zLD"]:
                sol = {"zLD":self.problem.value,"status":self.problem.status,"X":self.vp["X"].value,"lam":self.vp["lam"].value}
            
            values.append(self.problem.value)
            # if self.problem.value == sol["zLD"]:
            #     no_change_counter += 1
            #     if no_change_counter > 3:
            #         if betha == 2:
            #             betha = 1
            #             no_change_counter = 0
            #         else:
            #             print("no change in lambda limit reached in ", t,"-th step - converged")
            #             sol["converged"] = t
            #             break  
            # else:
            #     no_change_counter = 0
            
            
            # subgradient = prosta kapaciteta
            s = self.vp["cap"].value - np.sum(self.vp["X"].value,axis=1)
            if np.linalg.norm(s) == 0:
                raise("miracle: s == 0")
                return
            # izberemo korak (nalivno ali op)
            # 1
            # alpha = 1/(t+1)
            
            # 2
            # if UB == np.inf: 
            #     alpha = 1/(t+1)
            # else:
            #     alpha = 1/(t+1) * (UB - self.problem.value)/np.linalg.norm(s)**2  # TODO
            
            # 3
            
            if UB == np.inf: 
                alpha = 1/(t+1)
            else:
                alpha = betha * (UB - self.problem.value)/np.linalg.norm(s)**2  # TODO
            
            # print("lambda:",self.vp["lam"].value)
            # print("s:",s)
            # print("alpha:",alpha)
            ll = self.vp["lam"].value - s * alpha
            ll[ll < 0] = 0 # lambda ne mora biti negativna
            # if np.all(self.vp["lam"].value == ll):
            #     no_change_counter += 1
            #     if no_change_counter > 3:
            #         print("no change in lambda limit reached in ", t,"-th step - converged")
            #         sol["converged"] = t
            #         break
            # else:
            #     no_change_counter = 0
            self.vp["lam"].value = ll
            
        
        #hf.plot_solution_graph(graph,sol["X"])
        # print("status:", sol["status"])
        # print("optimal value zLD", sol["val"])
        # print("lambda:", sol["lam"])
        # print("optimal var x", sol["X"])
        sol["s"] = self.vp["cap"].value - np.sum(sol["X"],axis=1)
        sol["z"] = self.vp["c"].value.T @ np.sum(sol["X"],axis=1)
        sol["cap_ok"] = np.all(np.sum(sol["X"],axis=1) <= self.vp["cap"].value)
        sol["zLD_ceil"] = np.ceil(sol["zLD"])
        
        self.sol = sol
        
        return values
    
    # def solve(self,LB_,UB_,lam0 = None,MAX_ITER_LR=100):
    #     self.label_solved = Node.label_solved
    #     Node.label_solved += 1
        
    #     sol = {"zLD":-np.inf,"status":None,"X":None,"lam":None,"cap_ok":False}
        
    #     if lam0 is None:
    #         self.vp["lam"].value = np.zeros(self.vp["lam"].value.shape)
    #     else:
    #         self.vp["lam"].value = lam0
        
    #     q = 0 # initial number of iteration
    #     flag = 0
    #     betha = 2
    #     q_max = MAX_ITER_LR # max number of iteration is q_max.
    #     LB = LB_
    #     UB = UB_
    #     # UB = np.sum(vp["c"].value) * vp["H"].value.shape[1] * np.max(vp["H"].value) # to je vsi komoditiji grejo po vseh povezavah
    #     # LB = -np.inf  # initial upper bound and lower bound 
    #     eps = 10**(-11)
        
    #     UB_min = UB
    #     X_best = None
        
    #     # problem = cp.Problem(obj, constraints)
        
    #     # Node.label = 0
    #     # Node.INIT_NUM_STEPS = INIT_NUM_STEPS
    #     # n = Node(obj,constraints,vp)
    #     values = []
    #     while q <= q_max and betha > eps:
    #         self.problem.solve()
    #         X = self.vp["X"].value
    #         zLD = self.problem.value
    #         values.append(zLD)
    #         if np.all(np.sum(X,axis=1) <= self.vp["cap"].value):#X* is feasible:
    #             UB = self.vp["c"].value.T @ np.sum(X,axis=1) # z
    #             if UB < UB_min:
    #                 UB_min = UB
    #                 sol["X"] = X
                    
    #         if zLD < LB:
    #             flag = 3
    #         else:
    #             if zLD - LB < eps * max(1,LB):
    #                 flag = flag + 1
    #             if zLD > LB:
    #                 LB = zLD
                    
    #         if flag > 2:
    #             betha = betha/2
    #             flag = 0

    #         s = self.vp["cap"].value - np.sum(self.vp["X"].value,axis=1)
    #         alpha = betha * (UB - zLD)/np.linalg.norm(s)
    #         ll = self.vp["lam"].value - s * alpha
    #         ll[ll < 0] = 0 # lambda ne mora biti negativna
    #         self.vp["lam"].value = ll
    #         q = q + 1
    #     print(q,betha)
        
    #     sol["s"] = self.vp["cap"].value - np.sum(sol["X"],axis=1)
    #     sol["z"] = UB # self.vp["c"].value.T @ np.sum(sol["X"],axis=1)
    #     sol["cap_ok"] = np.all(np.sum(sol["X"],axis=1) <= self.vp["cap"].value)
    #     sol["zLD_ceil"] = LB # np.ceil(sol["zLD"])
        
    #     self.sol = sol
        
    #     return values
    
    
    
    


    
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
        
        r = np.argsort(-s)
        for a in r:
            k = np.random.randint(0,X.shape[1])
            val = int(X[a,k])
            branching = (a,k,val,val+1)
            if is_new_branching(branching):
                break
            
        
        constraints1 = self.problem.constraints + [self.vp["X"][a,k] <= val] # false
        constraints2 = self.problem.constraints + [self.vp["X"][a,k] >= (val +1)] # true        
        
        
        ch1 = Node(self.problem.objective,constraints1,self.vp,level=self.level+1,parent=self,branchingTF=branching + (False,))
        ch2 = Node(self.problem.objective,constraints2,self.vp,level=self.level+1,parent=self,branchingTF=branching + (True,))
        
        self.children = [ch1,ch2]
        
        return self.children
    
    def __str__(self, level=0):
        ret = "\t"*level+repr(self)+"\n"
        if self.children is not None:
            for child in self.children:
                ret += child.__str__(level+1)
        return ret
    
    def __repr__(self):
        if self.sol is None: return str(self.label) +": not available"
        else: stri =  str(self.label) + " / " + str(self.label_solved)
        
        if self.sol["status"] == "infeasible": return stri +": "+ self.sol["status"]
        stri = (stri +": "+ self.sol["status"]+" z: "+ str(self.sol["z"])+" zLD: "+ # str(self.sol["zLD"]) +
            "("+str(self.sol["zLD_ceil"])+")"+ " cap_ok = " + str(self.sol["cap_ok"]))
        # if "vejanje" in self.sol: stri += " " + repr(self.sol["vejanje"]) + " "
        if self.children is not None: stri += " " + repr(self.children[0].branchingTF[:4]) + " "
        # if "converged" in self.sol: stri += " zLD CONVERGED " + str(self.sol["converged"])
        return stri
