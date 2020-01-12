import numpy as np
from scipy.optimize import fsolve
from scipy.integrate import quad
import matplotlib.pyplot as plt
import math
import csv

Cpass = 1.043
Cpassd = 8.293
Cact = 1.596
Cactd = 0.6804
Cactdd = 0.2905
Cmyo = 10.1
Cshear = 0.258
Cmeta = 30
Ctoned = -5.222
Ctonedd = 10.11
D0 = 156.49 * 1e-6

shear = 0
h = np.exp(Ctoned)
# P = 100*133.333
Diameter = []
Pressure = []
f1 = []
f2 = []
Tes = []
vis = 0.7 * 1.0e-3
Q = (1 / 3) * 1e-09

sign = lambda x: math.copysign(1, x)

#Coefficients for metabolic response
d_tissue = 18.8*1e-6
M_0 = 8.28
c0 = 0.5
H_D = 0.4
H_T = 0.3
R_0 = 1.4e-3
R_1 = 0.891
zeta_H = 2.7
P_50 = 26.0*133.322
S_0 = 0.97
C_0 = 0.5
k_d = 2.0e-6
tau_d = 1
tau_a = 60
L_0 = 1.0e-2
xmp_La = 0.90e-2
y_L = 2.57e-2


def Tension2(D, P):
    
    meta = 0
    D_a = D+d_tissue
    q1 = (np.pi)*((D_a*D_a)-(D*D))*M_0
    alpha = ((H_T*R_0)/(4*k_d))*((D*(1-(R_1*S_0)))-(((1-H_D)*R_1*q1)/((np.pi)*c0*H_D*k_d)))
    beta = (D*H_T*R_0*R_1*q1)/(4*Q*c0*H_D*k_d)
    gamma = (k_d*(np.pi)*D)/(0.6*Q)
    add1 = C_0-alpha
    first_term = (alpha+beta*xmp_La)*L_0*(1-np.exp(-(y_L-xmp_La)/L_0))
    second_term = (L_0*L_0)*beta*(1-((np.exp(-(y_L-xmp_La)/L_0))*(1+((y_L-xmp_La)/L_0))))
    third_term = (add1*np.exp(-(gamma*xmp_La))/(gamma+(1/L_0)))*(1-np.exp(-(y_L-xmp_La)/(gamma+(1/L_0))))
    meta = first_term+second_term+third_term
    lam_D = D / D0
    shear = (32 * Q * vis) / ((np.pi) * math.pow(D, 3))
    Stone = (Cmyo * (P * D * 0.5)) - Ctonedd - (Cshear * shear)-(Cmeta*meta)
    Act   = 1 / (1 + np.exp(-Stone))
    p1    = (lam_D - 1) * Cpassd
    # print(p1)
    p2 = np.exp(p1)
    Tpass = Cpass * p2
    s1   = (lam_D - Cactd) / (Cactdd)
    s2   = math.pow(s1, 2)
    s3   = np.exp(-s2)
    Tact = Cact * s3
    Ttot = Tpass + (Act * Tact)
    return (Ttot) - (P * D * 0.5)

def O2_consumption(x):
    
    x0 = 0
    if(x<0.61):
        return 0.5
    elif(0.61<x<=1.2):
        D = 65.2*1e-6
        x0 = 0
        n_i = 1
    elif(1.2<x<=1.56):
        D = 14.8*1e-6
        n_i = 143
        x0 = 0.59*1e-2
    elif(1.56<x<=1.61):
        D = 6*1e-6
        n_i = 5515
        x0 = 0.95*1e-2
    elif(1.61<x<=1.97):
        D = 23.5*1e-6
        x0 = 1e-2
        n_i = 143
    elif(1.97<x<=2.56):
        D = 119.1*1e-6
        x0 = 1.95*1e-2
        n_i = 1
    else:
        return 2
    x = x*1e-2
    D_a = D +(2*d_tissue)
    q1 = (np.pi) * ((D_a * D_a) - (D * D)) * M_0*0.25*0.01
    alpha = ((H_T * R_0) / (4 * k_d)) * ((D * (1 - (R_1 * S_0))) - (((1 - H_D) * R_1 * q1) / ((np.pi) * c0 * H_D * k_d)))
    beta = (D * H_T * R_0 * R_1 * q1) / (4 * Q * c0 * H_D * k_d)
    gamma = (k_d * (np.pi) * D) / (0.6 * Q)
    C_x = alpha + (beta*x) + ((np.exp(gamma*(x0-x)))*(C_0-alpha-beta*x0))
    return C_x

def Saturation(x):
    
    x0 = 0
    if(x<=0.61):
        return 0.97
    elif(0.61<x<=1.2):
        D = 65.2*1e-6
        x0 = 0.61
        n_i = 1
        vel = 2.13*1e-2
    elif(1.56<x<=1.61):
        D = 6*1e-6
        vel = 0.05*1e-2
        x0 = 1.56
        n_i = 5515
    elif(1.2<x<=1.56):
        D = 14.8*1e-6
        vel = 0.29*1e-2
        x0 = 1.2
        n_i = 143
    elif(1.61<x<=1.97):
        return Saturation(1.61)
    elif(1.97<x<=2.56):
        return Saturation(1.61)
    else:
        return Saturation(1.61)
    x1 = (x)*1e-2
    D_a = D + (2*d_tissue)
    Q  =  5.33*1e-8
    q1 =  (np.pi) * ((D_a * D_a) - (D * D)) * M_0*0.25*0.01
    S_x = Saturation(x-0.1) - ((q1*n_i)/(Q*c0*H_D))*x1
    if(S_x<0):
        return 0.426
    return S_x

def Activation(D,P):
    
    meta = 0
    D = D*1e-6
    P = P*133
    D_a = D + d_tissue
    q1 = (np.pi) * ((D_a * D_a) - (D * D)) * M_0
    alpha = ((H_T * R_0) / (4 * k_d)) * ((D * (1 - (R_1 * S_0))) - (((1 - H_D) * R_1 * q1) / ((np.pi) * c0 * H_D * k_d)))
    beta = (D * H_T * R_0 * R_1 * q1) / (4 * Q * c0 * H_D * k_d)
    gamma = (k_d * (np.pi) * D) / (0.6 * Q)
    add1 = C_0 - alpha
    first_term = (alpha + beta * xmp_La) * L_0 * (1 - np.exp(-(y_L - xmp_La) / L_0))
    second_term = (L_0 * L_0) * beta * (1 - ((np.exp(-(y_L - xmp_La) / L_0)) * (1 + ((y_L - xmp_La) / L_0))))
    third_term = (add1 * np.exp(-(gamma * xmp_La)) / (gamma + (1 / L_0))) * (1 - np.exp(-(y_L - xmp_La) / (gamma + (1 / L_0))))
    meta = first_term + second_term + third_term
    shear = (32 * Q * vis) / ((np.pi) * math.pow(D, 3))
    Smyo = (Cmyo * (P * D * 0.5))
    Sshear = (Cshear*shear)
    Smeta = (Cmeta*meta)
    Stone = Smyo - Sshear - Smeta - Ctonedd
    Act = 1 / (1 + np.exp(-Stone))
    #print(Smeta)
    return Act


D1 = range(5, 200, 5)
# Pres = range(40,95,5)
D_1 = range(1, 200, 1)
Diam = 0.001
D2 = [i / 150 for i in D1]
X0 = 0.001
k = 1
flag = 0
Pres = 5
D111 = 0
while (Pres <= 198):
    
    D111 = Pres * 133.3333
    # X0 = Dekkers(Tension2,0,0.00005,D111)
    Diam = fsolve(Tension2, Diam, args=D111)
    k = Tension2(Diam, D111)
    #print(k)
    if (abs(k) > 2):
        X0 = Diam + 0.000001
        print(Pres)
        continue

    Pressure.append(Pres)
    Diameter.append(Diam * 1e06)
    Pres += 1

Test_Pressure = []
Test_Diameter = []

for d in csv.DictReader(open('C:/Users/ASUS/Desktop/Photos/My_Work/Saturation.csv')):
    Test_Pressure.append(float(d['Distance']))
    Test_Diameter.append(float(d['Saturation']))


Stimuli = [Activation(ki,li) for ki,li in zip(Diameter,Pressure)]
x1 = np.arange(0,4,0.1).tolist()
x11 = list(map(lambda x:round(x,2),x1))
print(Saturation(1.3))
oxygen_st = [O2_consumption(li) for li in x1]
Sat = [Saturation(mi) for mi in x11]
plt.xlim(0,4)
plt.ylim(0,1)
print('Pressure = ', x11)
print('Diameter = ', Sat)
plt.xlabel('Pressure')
plt.ylabel('Activation')
plt.plot(Test_Pressure, Test_Diameter)
plt.plot(x11,Sat,'r')
plt.plot()
plt.show()


