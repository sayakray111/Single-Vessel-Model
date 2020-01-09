import matplotlib.pyplot as plt 
import numpy as np
from decimal import Decimal
X = []
Cp = []
Cx = []
Cy = []
Cz = []
C = []
#f = open("/hpc/sray036/VirtualPregnancy/repro-examples/uteroplacental/darcy-models/02-darcy-ell-full/output/StaticDarcy.part0.exnode",'r')
f = open("D:/Projects/repro-examples/uteroplacental/darcy-models/00-darcy-testcase/output/STATIC_SOLUTION.part0.exnode")
currentline  = 1
ST = f.readlines()

sum1 = 0.0
for n in range(1,7000):
    #x_val = float(ST[9*n+6].strip())
    x_val = round(float(ST[14*n+6].strip()),2)
    #y_val = float(ST[9*n+7].strip())
    y_val = round(float(ST[14*n+7].strip()),2)
    #z_val = float(ST[9*n+8].strip())
    z_val = round(float(ST[14*n+8].strip()),2)
    #print(x_val)
    if(y_val == 0.8 and z_val == 0.6):
        C.append(round((float(ST[14*n+11].strip())),2))
        X.append(x_val)
        #C1.append(np.exp(x_val+1))
        Cp.append(x_val*y_val*z_val)
        Cz.append(-x_val*y_val)
        Cy.append(-x_val*z_val)
        Cx.append(-z_val*y_val)

for k in range(0,len(Cz)):
    error = (C[k]-Cz[k])
    sum1 = sum1+error**2
    avg = sum1/len(C)
    avg = np.sqrt(avg)
print('Average error is:', avg)
plt.plot(X,C)
plt.plot(X,Cz)
plt.show()
f.close()