import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
import math
import csv
Cpass = 1.043
Cpassd = 8.293
Cact = 1.596
Cactd = 0.6805
Cactdd = 0.2905
Cmyo = 11.384
Cshear = 0.1
Cmeta = 30
Ctoned = -5.64066865
Ctonedd = 10.11
D0 = 156.49*1e-6

h = np.exp(Ctoned)
#P = 100*133.333
Diameter = []
Pressure = []
f1  = []
f2  = []
Tes = []
vis = 2.11*1e-3
Q = 1e-09
sign = lambda x: math.copysign(1, x)


def Tension2(D,P):
    lam_D = D/D0
    shear = (32 * Q * vis) / ((np.pi) * math.pow(D, 3))
    #shear = 0
    Stone = (Cmyo*(P*D*0.5))+Ctoned-(Cshear*shear)
    Act = 1/(1+np.exp(-Stone))
    p1 = (lam_D-1)*Cpassd
    #print(p1)
    p2 = np.exp(p1)
    Tpass = Cpass * p2
    s1 = (lam_D-Cactd)/(Cactdd)
    s2 = math.pow(s1,2)
    s3 = np.exp(-s2)
    Tact = Cact * s3
    Ttot = Tpass+(Tact*Act)
    #print('Activation = ',Act,'Passive Tension = ',Tpass,'Active Tension = ',(Tact*Act))
    #print(Ttot-(P*D*0.5))
    return (Ttot)-(P*D*0.5)



def Activation(P,D):
    T = P*D*0.5
    Stim = (Cmyo*T)+Ctoned
    return 1/(1+np.exp(-Stim))

def Dekkers(f,a,b,d):
    fa = f(a,d)
    fb = f(b,d)
    k = 1
    if(fa*fb>0):
        #print(fa,'a value',fb)
        return 1
    c = a
    fc = fa
    while(True):
        if(sign(fb)==sign(fc)):
            c = a
            fc = fa
        if(abs(fc)<abs(fb)):
            a,fa = b,fb
            b,fb = c,fc
            c,fc = a,fa
            
        m = (b+c)/2
        if(abs(m-b)<=np.spacing(abs(b))):
            return 
        p = (b-a)*fb
        if(p>=0):
            q = fa - fb
        else:
            q = fb - fa
            p = -p
        a = b
        fa = fb
        k = k+1
        
        if(p<=np.spacing(q)):
            b = b + sign(c-b) * np.spacing(b)
            fb = f(b,d)
            root = b
            return root
        elif(p<=(m-b)*q):
            b = b + p/q
            fb = f(b,d)
            root = b
            return root
        else:
            b = m
            fb = f(m,d)
            root = b
            return root

D1 = range(5,200,5)
#Pres = range(40,95,5)
D_1 = range(1,200,1)
Diam = 0.001
Act = []
D2 = [i/150 for i in D1]
X0=0.001
k = 1
flag = 0
Pres = 5
D111 = 0
while(Pres<=150):
    D111 = Pres*133
    #X0 = Dekkers(Tension2,0,0.00005,D111)
    Diam = fsolve(Tension2,Diam,args=D111)
    k = Tension2(Diam,D111)
    k1 = Activation(D111,Diam)
    Act.append(k1)
    if(abs(k)>1e-10):
        X0 = Diam+1e-6
        continue

    Pressure.append(Pres)
    Diameter.append(Diam*1e06)
    Pres+=1
    


Test_Pressure = []
Test_Diameter = []

for d in csv.DictReader(open('Carlson(2008)_60l.csv')):
    Test_Pressure.append(float(d['Pressure']))
    Test_Diameter.append(float(d['Diameter']))

print('Pressure = ', Test_Pressure)
print('Diameter = ', Test_Diameter)
T1 = Test_Pressure[3]*Test_Diameter[3]*133*1e-6*0.5



plt.xlim(0,120)
plt.ylim(0,180)

plt.xlabel('Pressure')
plt.ylabel('Diameter')
plt.plot(Test_Pressure,Test_Diameter)
plt.plot(Pressure,Diameter)
#plt.xlim(0,120)
#plt.ylim(0,1)

#plt.xlabel('Pressure')
#plt.ylabel('Activation')
#plt.plot(Pressure,Act)
plt.show()

        
        