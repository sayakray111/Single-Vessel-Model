import numpy as np
from scipy.optimize import root_scalar,fsolve
import matplotlib.pyplot as plt
import math
import csv
import scipy.io
import scipy.integrate as integrate

mat = scipy.io.loadmat('DAnometa2vals.mat')
Pvc = float(mat.get('Pvc') / 1333)
Pac = float(mat.get('Pac') / 1333)
# calculate the value of Q to be used for the evaluation of shear stress
Q_tot = 0
gradP_tot = 0
Resistance_total = 0
Resistance_c = 0
Q_la = 0
Q_sa = 0
Q_c = 0
Q_sv = 0
Q_lv = 0
Q_la = Q_tot
Q_sa = Q_tot / 143
Q_lv = Q_tot
Q_sv = Q_tot / 143
Q_c = Q_tot / 5515

# calculate the value of Diameter of segments....
D_c = 6 * 1e-6
D_la = 0
D_lac = 65.2 * 1e-6
D_sac = 14.8 * 1e-6
D_sa = 0
D_lv = 23.5 * 1e-6
D_sv = 119.1 * 1e-6

# The length of the various segments.

l_c = 0.05 * 0.01
l_v = 0.61 * 0.01
l_a = 0.61 * 0.01
l_la = 0.59 * 0.01
l_lv = 0.59 * 0.01
l_sv = 0.36 * 0.01
l_sa = 0.36 * 0.01
# the pressure drop in the segments.
P_la = 0
P_lv = 0
P_sa = 0
P_sv = 0
P_c = 0
# Numbers of various vessels.
n_c = 5515
n_la = 1
n_a = 1
n_v = 1
n_sa = 143
n_lv = 1
n_sv = 143
# volume of the segments...
vol_la = 0
vol_sa = 0
vol_c = 0
vol_sv = 0
vol_lv = 0
vol_tot = 0
shear = 0

# P = 100*133.333
Diameter = []
Pressure = []
f1 = []
f2 = []
Tes = []
vis = 7 * 1e-4

sign = lambda x: math.copysign(1, x)
N1 = 0


def Tension2(D, *params):
    
    Pres,Cpass,Cpassd,Cact,Cactd,Cactdd,Cmyo,Cshear,Cmeta,Ctonedd,D0,Diam_2,flag,N1,vis_1,xD = params
    # Determine which diameter is which....
    if(flag==1):
        Diam_la = D
        Diam_sa = Diam_2
    else:
        Diam_la = Diam_2
        Diam_sa = D
    #Calculate the Discharge for the given diameters and their initial guesses...
    gradP_tot = (Pres - 14) * 133.33
    Resistance_la = (128 * vis_la * l_la) / (np.pi * math.pow(Diam_la, 4) * n_la)
    Resistance_sa = (128 * vis_sa * l_sa) / (np.pi * math.pow(Diam_sa, 4) * n_sa)
    Resistance_c  = 2.81*1e13
    Resistance_lv = 0.86*1e13
    Resistance_sv = 0.28*1e13
    Resistance_a =  (128 * vis_a * l_a) / (np.pi * math.pow(Diam_a, 4) * n_a)
    Resistance_v =  (128 * vis_v * l_v) / (np.pi * math.pow(Diam_v, 4) * n_v)
    Resistance_total = Resistance_sa + Resistance_la + Resistance_c + Resistance_lv + Resistance_sv
    Resistance_total = Resistance_total + Resistance_a
    Q_tot = gradP_tot / Resistance_total
    #-----------------------------------------------------
    Q = Q_tot/N1
    if(flag==1):
        shear = (32 * Q * vis_1)/((np.pi) * math.pow(Diam_la, 3))
        P1 = (Pres * 133.333)-(Q_tot*Resistance_a)
        P2 = (Pres * 133.333)- (Q_tot * Resistance_la)-(Q_tot*Resistance_a)
        Pmid = (P1+P2)*0.5
    else:
        shear = (32 * Q * vis_1)/((np.pi) * math.pow(Diam_sa, 3))
        P2 = (Pres * 133.333)- (Q_tot * Resistance_la)-(Q_tot*Resistance_sa)-(Q_tot*Resistance_a)
        P1 = (Pres * 133.333)- (Q_tot * Resistance_la)-(Q_tot*Resistance_a)
        Pmid = (P1+P2)*0.5
    cons = (Diam_la,Diam_sa,Q_tot)
    meta = SCR(xD,*cons)
    lam_D = D / D0
    Stone = (Cmyo * (Pmid * D * 0.5)) + Ctonedd - (Cshear * shear) - (Cmeta * meta)
    Act = 1 / (1 + np.exp(-Stone))
    p1 = (lam_D - 1) * Cpassd
    p2 = np.exp(p1)
    Tpass = Cpass * p2
    s1 = (lam_D - Cactd) / (Cactdd)
    s2 = math.pow(s1, 2)
    s3 = np.exp(-s2)
    Tact = Cact * s3
    Ttot = Tpass + (Act * Tact)
    return (Ttot) - (Pmid * D * 0.5)



def Saturation(x,*params1):
    D1,D2,Q_tot = params1
    Pw2 = (D1,D2,Q_tot)
    if (0 < x <= 0.61):
        return 0.97
    elif (0.61 < x <= 1.2):
        D = D1
        S_i = 0.97
        Q = Q_tot
        X_i = 0.61 * 1e-2
    elif (1.2 < x <= 1.56):
        D = D2
        S_i = Saturation(1.2,*Pw2)
        Q = Q_tot / 143
        X_i = 1.2 * 1e-2
    elif (1.56 < x <= 1.61):
        D = 6 * 1e-6
        S_i = Saturation(1.56,*Pw2)
        Q = Q_tot / 5515
        X_i = 1.56 * 1e-2
    elif (1.61 < x <= 1.97):
        return Saturation(1.61,*Pw2)
    else:
        return Saturation(1.97,*Pw2)
    q1 = 0.25 * np.pi * M_0 * (((D + (37.6 * 1e-6)) * (D + (37.6 * 1e-6))) - (D * D))
    x = x * 1e-2
    A1 = S_i - (q1 / (Q * c0 * H_D)) * (x - X_i)
    return S_i - (q1 / (Q * c0 * H_D)) * (x - X_i)


def Consumption(x,*params2):
    
    D1,D2,Q_tot  = params2
    Pw = (D1,D2,Q_tot)
    X = x * 1e-2
    if (0 < x <= 0.61):
        C_i = 0.5 * 1e-3
        return C_i
    elif (0.61 < x <= 1.2):
        D = D1
        C_i = 0.5 * 1e-3
        S_i = 0.97
        Q = Q_tot
        X_i = 0.61 * 1e-2
    elif (1.2 < x <= 1.56):
        D = D2
        C_i = Consumption(1.2,*Pw)
        S_i = Saturation(1.2,*Pw)
        Q = Q_tot / 143
        X_i = 1.2 * 1e-2
    elif (1.56 < x <= 1.61):
        D = 6 * 1e-6
        C_i = Consumption(1.56,*Pw)
        S_i = Saturation(1.56,*Pw)
        Q = Q_tot / 5515
        X_i = 1.56 * 1e-2
    elif (1.61 < x <= 1.97):
        D = 23.5 * 1e-6
        C_i = Consumption(1.61,*Pw)
        S_i = Saturation(1.61,*Pw)
        Q = Q_tot / 143
        X_i = 1.61 * 1e-2
        gamma1 = (k_d * (np.pi) * D) / (0.6 * Q)
        alpha1 = ((H_T * R_0) / (4 * k_d)) * ((D * (1 - (R_1 * S_i))))
        A1 = (C_i - alpha1) / (np.exp(-gamma1 * X_i))
        C = A1 * np.exp(-gamma1 * X) + alpha1
        return C
    else:
        D = 119.1 * 1e-6
        C_i = Consumption(1.97,*Pw)
        S_i = Saturation(1.97,*Pw)
        Q = Q_tot
        X_i = 1.97 * 1e-2
        gamma1 = (k_d * (np.pi) * D) / (0.6 * Q)
        alpha1 = ((H_T * R_0) / (4 * k_d)) * ((D * (1 - (R_1 * S_i))))
        A1 = (C_i - alpha1) / (np.exp(-gamma1 * X_i))
        C = A1 * np.exp(-gamma1 * X) + alpha1
        return C
    q1 = 0.25 * np.pi * M_0 * (((D + 37.6 * 1e-6) * (D + 37.6 * 1e-6)) - (D * D))
    alpha = ((H_T * R_0) / (4 * k_d)) * (
                (D * (1 - (R_1 * S_i))) - (((1 - H_D) * R_1 * q1) / ((np.pi) * c0 * H_D * k_d)))
    beta = (D * H_T * R_0 * R_1 * q1) / (4 * Q * c0 * H_D * k_d)
    gamma = (k_d * (np.pi) * D) / (0.6 * Q)
    C = alpha + (beta * (X - X_i)) + np.exp(gamma * (X_i - X)) * (C_i - alpha)
    return C

def Consumption2(x,*params3):
    x11,D1,D2,Q_tot = params3
    params2 = (D1,D2,Q_tot)
    C = np.exp(-(x-x11)/(L_0))*Consumption(x*100,*params2)
    return C
def SCR(x,*params4):
    x1 = x*1e-2
    params41 = (x1,*params4)
    SCR1,error = integrate.quad(Consumption2,x1,xend_lv,args=params41)
    #print('error is = ',error)
    return SCR1
# Viscosity values....
vis_a = 0.0226 * 0.1
vis_la = 0.0211 * 0.1
vis_sa = 0.0355 * 0.1
vis_c = 0.0905 * 0.1
vis_sv = 0.0256 * 0.1
vis_lv = 0.0234* 0.1
vis_v = 0.0250 * 0.1
# Diameters of the various segments.....
Resistance_a = 1.88*1e13
G = (128*vis_a*l_a)/(np.pi*Resistance_a*n_a)
Diam_a = math.pow(G,0.25)
print('Diameter of the main artery is = ',Diam_a*1e06)
Resistance_v = 0.19*1e13
G = (128*vis_v*l_v)/(np.pi*Resistance_v*n_v)
Diam_v = math.pow(G,0.25)
print('Diameter of the main vein is = ',Diam_v*1e06)
Diam_c = 6*1e-6
Diam_lv = 119.1*1e-6
Diam_sv = 23.5*1e-6
Diam_sa = 0.0001
Diam_la = 0.001
d_t = 18.8 * 1e-6
Diam_lac = 65.2*1e-6
Diam_sac = 14.8*1e-6
Diam_la1 = Diam_la + (2 * d_t)
Diam_lac1 = Diam_lac + (2 * d_t)
Diam_sa1 = Diam_sa + (2 * d_t)
Diam_sac1 = Diam_sac + (2 * d_t)
Diam_c1 = Diam_c + (2 * d_t)
Diam_sv1 = Diam_sv + (2 * d_t)
Diam_lv1 = Diam_lv + (2 * d_t)
#######
perfuse = 1
perfusion = []
perfuse_100 = 0
Pressure_sa = []
Diameter_sa = []
Diameter_la = []
Pressure_la = []
Pres = 51
D111 = 0.0
perfuse = 0.0
metac = 0.0
Pressure_in = []
# Coefficients for metabolic response... small arteriole. large arteriole...
d_tissue = 18.8 * 1e-6
M_0 = (8.28 / 60) * 0.01
c0 = 0.5
H_D = 0.4
H_T = 0.3
R_0 = 1.4e-3
R_1 = 0.891
zeta_H = 2.7
P_50 = 26.0 * 133.322
S_0 = 0.97
C_0_la = 0.5 * 1e-3
k_d = 2.0 * 1e-6
tau_d = 1
tau_a = 60
L_0 = 1.0e-2
x0_la = 0.61e-2
x0_sa = 1.2e-2
x0_c = 1.56e-2
x0_sv = 1.61e-2
x0_lv = 1.97e-2
x0_v = 2.56e-2
xmp_la = 0.90e-2
xmp_sa = 1.38e-2
xend_la = 1.2e-2
xend_sa = 1.56e-2
xend_c = 1.61e-2
xend_sv = 1.97e-2
xend_lv = 2.56e-2
d_t = 18.8 * 1e-6
# Parameters for large arteriole....
Cpass_la = 1.043
Cpassd_la = 8.293
Cact_la = 1.596
Cactd_la = 0.6804
Cactdd_la = 0.2905
Cmyo_la = 10.1
Cshear_la = 0.258
Cmeta_la = 3 * 1e06
Ctoned_la = -5.22
Ctonedd_la = 10.11
D0_la = 156.49 * 1e-6
 # small arteriole...
Cpass_sa = 0.2599
Cpassd_sa = 11.467
Cact_sa = 0.274193
Cactd_sa = 0.750
Cactdd_sa = 0.384
Cmyo_sa = 35.9
Cshear_sa = 0.258
Ctoned_sa = -15.3
Cmeta_sa = 3 * 1e06
D0_sa = 38.99 * 1e-6
Ctonedd_sa = 10.66

Activation_la = []
Pvc = 11.91
#Pvc = 45 ####Only for the shear case
Pres = 20 ##Only for the myogenic case

while (Pres <= 200):
    # large arteriole...
    params = (Pres,Cpass_la,Cpassd_la,Cact_la,Cactd_la,Cactdd_la,Cmyo_la,Cshear_la,Cmeta_la,Ctonedd_la,D0_la,Diam_sa,1,n_la,vis_la,0.90)
    Diam_la = fsolve(Tension2,Diam_la,args=params)
    k = Tension2(Diam_la, *params)
    if(abs(k)>1e-2 or Diam_la<1e-6):
        n  = 1
        while(True):
            Diam_la = (Diameter_la[-n]*1e-6)+(n*1e-8)

            Diam_la = fsolve(Tension2,Diam_la,args=params)
            print('Hitting a hole at a diameter (large) = ', Diam_la, ' Pressure  = ', Pres)
            k = Tension2(Diam_la, *params)
            print('error is = ',k)
            if(Diam_la>3.5e-5 and abs(k)<0.1):
                break
            n+=1
    d_t = 18.8 * 1e-6
    #small arteriole...
    params = (Pres,Cpass_sa,Cpassd_sa,Cact_sa,Cactd_sa,Cactdd_sa,Cmyo_sa,Cshear_sa,Cmeta_sa,Ctonedd_sa,D0_sa,Diam_la,2,n_sa,vis_sa,1.38)
    Diam_sa = fsolve(Tension2,Diam_la,args=params)
    k = Tension2(Diam_sa, *params)
    if(Diam_sa<1e-6 or abs(k)>1e-2):
        n = 1
        while(True):
            Diam_sa = (Diameter_sa[-1]*1e-6)+(n*1e-7)
            Diam_sa = fsolve(Tension2,Diam_sa,args=params)
            print('Hitting a hole at a diameter (small) = ', Diam_sa, ' Pressure = ', Pres)
            k = Tension2(Diam_sa, *params)
            print('error is small = ',k)
            if(Diam_sa<3e-5 and abs(k)<0.1):
                break
            n+=1
    Pressure_in.append(Pres)
    Pressure_la.append(P_la / 133.33)
    Diameter_la.append(Diam_la * 1e06)
    Pressure_sa.append(P_sa / 133.33)
    Diameter_sa.append(Diam_sa * 1e06)

    Pres+=1

k = 0
while(k<len(Diameter_la)):
    gradP_tot = (Pressure_in[k] - 14) * 133.33
    Resistance_la = (128 * vis_la * l_la) / (np.pi * math.pow(Diameter_la[k]*1e-6, 4) * n_la)
    Resistance_sa = (128 * vis_sa * l_sa) / (np.pi * math.pow(Diameter_sa[k]*1e-6, 4) * n_sa)
    Resistance_c = 2.81 * 1e13
    Resistance_lv = 0.86 * 1e13
    Resistance_sv = 0.28 * 1e13
    Resistance_a = (128 * vis_a * l_a) / (np.pi * math.pow(Diam_a, 4) * n_a)
    Resistance_v = (128 * vis_v * l_v) / (np.pi * math.pow(Diam_v, 4) * n_v)
    Resistance_total = Resistance_sa + Resistance_la + Resistance_c + Resistance_lv + Resistance_sv
    Resistance_total = Resistance_total + Resistance_a
    print(Resistance_total)
    Q_tot = gradP_tot / Resistance_total
    Diam_a1 = Diam_a + (2 * d_t)
    Diam_c1 = Diam_c + (2 * d_t)
    Diam_la1 = Diam_la + (2 * d_t)
    Diam_sa1 = Diam_sa + (2 * d_t)
    vol_a = np.pi * 0.25 * (Diam_a1 * Diam_a1) * l_a * n_a
    vol_la = np.pi * 0.25 * (Diam_la1 * Diam_la1) * l_la * n_la
    vol_sa = np.pi * 0.25 * (Diam_sa1 * Diam_sa1) * l_sa * n_sa
    vol_c = np.pi * 0.25 * (Diam_c1 * Diam_c1) * l_c * n_c
    vol_lv = np.pi * 0.25 * (Diam_lv * Diam_lv) * l_lv * n_lv
    vol_sv = np.pi * 0.25 * (Diam_sv * Diam_sv) * l_sv * n_sv
    vol_v = np.pi * 0.25 * (Diam_v * Diam_v) * l_v * n_v
    vol_tot = vol_la + vol_sa + vol_c + vol_lv + vol_sv + vol_a + vol_v
    perfuse = Q_tot / (vol_tot)
    perfusion.append(perfuse)

    if (Pressure_in[k] == 100):
        perfuse_100 = perfuse
        print('Perfusion=', perfuse_100)

    k+=1


#perfuse_100 = (70.9)/(6000*6000)
Diam_lac = 65.2 * 1e-6
vis_la = 0.00211
l_la = 0.59 * 1e-2
n_la = 1
Resistance_lac = (128 * vis_c * l_c) / (np.pi * math.pow(Diam_c, 4) * n_c)
print('Value got = ', Resistance_lac / (1e13))
print('Value in paper = ', (7.5))
l1 = np.arange(0.1, 3, 0.01).tolist()
perfusion_norm = [(k / perfuse_100) for k in perfusion]
#Sat = [Saturation(l) for l in l1]
#Conc = [Consumption(l) * 1000 for l in l1]
Test_Pressure = []
Test_Pressure1 = []
Test_Diameter = []
Test_perfusion = []
#gradP_100 = ((70.9 / 6000) * (Resistance_total100) * vol_100) / 133
#print(Diameter_sa)
#print(S3)
for d in csv.DictReader(open('Carlson_2008(myo+shear+meta).csv')):
    Test_Pressure.append(float(d['Pressure']))
    Test_Diameter.append(float(d['Diameter']))
for d in csv.DictReader(open('perfusion(myo+shear+meta).csv')):
    Test_Pressure1.append(float(d['Pressure']))
    Test_perfusion.append(float(d['Perfusion']))
# print('Pressure = ', Test_Pressure)
# print('Diameter = ', Test_Diameter)
plt.figure(figsize=(10,9))
plt.subplot(211)
plt.xlim(0, 200)
plt.ylim(0,180)
plt.plot(Test_Pressure, Test_Diameter, 'b')
plt.plot(Pressure_in,Diameter_la,'r')
plt.subplot(212)
plt.xlim(0, 200)
plt.ylim(0,3)
plt.plot(Test_Pressure1, Test_perfusion, 'b')
plt.plot(Pressure_in,perfusion_norm,'r')
plt.show()