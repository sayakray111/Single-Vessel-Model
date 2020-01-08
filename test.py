
import matplotlib.pyplot as plt 
import numpy as np
from decimal import Decimal
X = []
C = []
C1 = []
f = open("/hpc/sray036/opencmiss/02-darcy-ell-full/output/STATIC_SOLUTION.part0.exnode",'r')
#f = open("/hpc/sray036/VirtualPregnancy/repro-examples/uteroplacental/advection-diffusion-models/01-advdiff-cuboid/output/iron/My_Work.part0.exnode",'r')
ST = f.readlines()
print(len(ST))
for n in range(1,5000):

    x_val = float(ST[14*n+6].strip())
    y_val = float(ST[14*n+7].strip())

    z_val = float(ST[14*n+8].strip())

    if(y_val == 0.5 and z_val == 0.5):
        #C.append(float(ST[9*n+9].strip()))
        C.append(float(ST[14*n+12].strip()))

        X.append(x_val)
        #C1.append(np.exp(x_val+1))
        C1.append(x_val*0.25)

#print(X)
plt.plot(X,C)
plt.plot(X,C1)
plt.show()
f.close()
