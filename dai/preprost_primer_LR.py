from matplotlib import pyplot as plt
import numpy as np

c = np.array([10,5])
cap = np.array([2,1])
Q = [0,1,2]

# lam = np.zeros(2)
lam = np.array([0,10])
lams = []
zs = []
alphas = []
for i in range(100):
    print("lam: ",lam)
    a = c[0] +lam[0] - c[1] - lam[1]
    b = 2*(c[1]+lam[1]) - lam[0]*cap[0] - lam[1]*cap[1]
    
    lams.append(lam)
    
    x11 = min(Q, key=lambda x: (not x == 1) if abs(a) < 10e-7 else a * x)
    z = a*x11 + b
    zs.append(z)
    plt.scatter(lam[1],z)
    # plt.plot(Q, [a*x+b for x in Q])

    print("x11: ",x11)
    alpha = 1/(i+1)
    alphas.append(alpha)
    print(alpha)
    s = cap - np.array([x11,2-x11])
    print("s: ",s)
    lam = lam - alpha * s
    lam[lam<0] = 0
plt.title("približek za zLD v odvisnosti od lambde")
plt.xlabel("lambda")
plt.ylabel("g(lambda)")
plt.show()
# plt.title("conv(h) pri fiksnem lambda")
# plt.xlabel("x_1,1")
# plt.ylabel("conv(h)")
# plt.show()
plt.plot([lam[1] for lam in lams])
plt.xlabel("i. korak")
plt.ylabel("lambda_2")
plt.title("lambda")
plt.show()
plt.plot(zs) 
plt.title("približki zLD")
plt.xlabel("i. korak")
plt.ylabel("g(lamda^i)")
plt.show()
plt.plot(alphas) 
plt.title("velikost koraka")
plt.ylabel("alpha")
plt.xlabel("i. korak")
plt.show()
