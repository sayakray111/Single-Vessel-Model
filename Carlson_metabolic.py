import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
import math
import csv
import scipy.io
import scipy.integrate as integrate

def compartment_resistance(viscosity,length,diameter,number_in_generation):
	resistance  = (128. * viscosity* length) / (np.pi * math.pow(diameter, 4.) * number_in_generation)
	return resistance

# calculate the value of Q to be used for the evaluation of shear stress
Q_tot = 0
gradP_tot = 0
Resistance_total = 0
# The length of the various segments converted to metre from centimetre
# Adapted from the paper by arciero et.al 
l_c = 0.05 * 0.01
l_v = 0.61 * 0.01
l_a = 0.61 * 0.01
l_la = 0.59 * 0.01
l_lv = 0.59 * 0.01
l_sv = 0.36 * 0.01
l_sa = 0.36 * 0.01
# the pressure drop in the segments.
P1 = 0.0
P2 = 0.0
# Numbers of various vessels adapted from the paper....
n_c = 5515
n_la = 1
n_a = 1
n_v = 1
n_sa = 143
n_lv = 1
n_sv = 143
# volume of the segments...
vol_la = 0.0
vol_sa = 0.0
vol_c = 0.0
vol_sv = 0.0
vol_lv = 0.0
vol_tot = 0.0
shear = 0.0
meta = 0.0


# This function is input into fsolve to calculate the Optimum Diameter for a particular pressure....
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
    gradP_tot = (Pres - 14) * 133.33 # Pressure converted from mmHg to N/m^2 by multiplying with 133.33 End vein pressure assumed to be 14mmHg
    # Calculating the resistances here ....
    Resistance_la = compartment_resistance(vis_la,l_la,Diam_la,n_la) 
    Resistance_sa = compartment_resistance(vis_sa,l_sa,Diam_sa,n_sa) 
    Resistance_c  = compartment_resistance(vis_c,l_c,Diam_c,n_c) 
    Resistance_lv = compartment_resistance(vis_lv,l_lv,Diam_lv,n_lv) 
    Resistance_sv = compartment_resistance(vis_sv,l_sv,Diam_sv,n_sv) 
    Resistance_a =  compartment_resistance(vis_a,l_a,Diam_a,n_a) 
    Resistance_v =  compartment_resistance(vis_v,l_v,Diam_v,n_v) 
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
    # Calculating the consumption of ATP...
    cons = (Diam_la,Diam_sa,Q_tot)
    # Calculating the metabolic constant (SCR) ...
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


# This function returns the saturation of Oxygen as a function of length....
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

# This calculates the consumption of ATP as a function of x(distance along the length of the vessel)
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
# Function inside integral
def Consumption2(x,*params3):
    x11,D1,D2,Q_tot = params3
    params2 = (D1,D2,Q_tot)
    C = np.exp(-(x-x11)/(L_0))*Consumption(x*100,*params2)
    return C
# This function calculates the SCR value for each segment...
def SCR(x,*params4):
    x1 = x*1e-2
    params41 = (x1,*params4)
    SCR1,error = integrate.quad(Consumption2,x1,xend_lv,args=params41)
    #print('error is = ',error)
    return SCR1
# Viscosity values....
vis_a = 0.0226 * 0.1 #viscosity of artery (Pa.s), directly from Arciero et al
vis_la = 0.0211 * 0.1 #viscosity of large arteriole (Pa.s), directly from Arciero et al
vis_sa = 0.0355 * 0.1 #viscosity of small arteriole(Pa.s), directly from Arciero et al
vis_c = 0.0905 * 0.1 #viscosity of capillary (Pa.s), directly from Arciero et al
vis_sv = 0.0256 * 0.1 #viscosity of small venule (Pa.s), directly from Arciero et al
vis_lv = 0.0234* 0.1 #viscosity of large venule (Pa.s), directly from Arciero et al
vis_v = 0.0250 * 0.1 #viscosity of vein (Pa.s), directly from Arciero et al
# Diameters of the artery.....This is done to calculate the diameter of the artery and the vein from their resistances...
Resistance_a = 1.88*1e11*133
G = (128*vis_a*l_a)/(np.pi*Resistance_a*n_a)
Diam_a = math.pow(G,0.25)
print('Diameter of the main artery is = ',Diam_a*1e06)
# Diameters of the vein....This is done to calculate the diameter of the artery and the vein from their resistances...
Resistance_v = 0.19*1e11*133
G = (128*vis_v*l_v)/(np.pi*Resistance_v*n_v)
Diam_v = math.pow(G,0.25)
print('Diameter of the main vein is = ',Diam_v*1e06)
Diam_c = 6*1e-6 # Diameter of the capillary adapted from the paper 
Diam_lv = 119.1*1e-6 # Diameter of the large venule adapted from the paper
Diam_sv = 23.5*1e-6  # Diameter of the small venule adapted from the paper
Diam_sa = 0.0001 # Diameter of the small arteriole which is given as initial guess to the fsolve routine...
Diam_la = 0.001 # Diameter of the large arteriole which is given as initial guess to the fsolve routine...
d_t = 18.8 * 1e-6 # layer of tissue around the vessels...
# Control Diameters of the large arteriole adapted from the paper 
Diam_lac = 65.2*1e-6
# Control Diameters of the small arteriole adapted from the paper.
Diam_sac = 14.8*1e-6
# Here the diameters of the thicknesses of the tissue layers around the arteries are added to the original diameters..
# these values are never used except in the perfusion calculations.
Diam_la1 = Diam_la + (2 * d_t)
Diam_lac1 = Diam_lac + (2 * d_t)
Diam_sa1 = Diam_sa + (2 * d_t)
Diam_sac1 = Diam_sac + (2 * d_t)
Diam_c1 = Diam_c + (2 * d_t)
Diam_sv1 = Diam_sv + (2 * d_t)
Diam_lv1 = Diam_lv + (2 * d_t)
####### Initialise certain list and values 
perfuse = 1.0
perfusion = [] # List to store the values of the perfusion...
Diameter_sa = [] # list to store the solved diameters of the small arteriole...
Diameter_la = [] # list to store the solved diameters of the large arteriole...
Pressure_la = [] # list to store the input pressures...

D111 = 0.0
perfuse = 0.0
metac = 0.0
Pressure_in = []
# Coefficients for metabolic response... small arteriole. large arteriole... adapted from the paper.
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
# Coefficients for large arteriole.... adapted from the paper.
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
 # Coefficients for small arteriole... adapted from the paper.
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
 # Assumed to be 14 mmHg.

Pres = 20 ## For the passive myogenic shear and metabolic case...
while (Pres <= 200):
    # List of parameters of the large arteriole...
    params = (Pres,Cpass_la,Cpassd_la,Cact_la,Cactd_la,Cactdd_la,Cmyo_la,Cshear_la,Cmeta_la,Ctonedd_la,D0_la,Diam_sa,1,n_la,vis_la,0.90)
    # Solve for the diameter of the large arteriole...
    Diam_la = fsolve(Tension2,Diam_la,args=params)
    k = Tension2(Diam_la, *params)
    # Calculate the optimum diameter to reduce the number of fluctuations...
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
    
    # List of parameters for small arteriole...
    params = (Pres,Cpass_sa,Cpassd_sa,Cact_sa,Cactd_sa,Cactdd_sa,Cmyo_sa,Cshear_sa,Cmeta_sa,Ctonedd_sa,D0_sa,Diam_la,2,n_sa,vis_sa,1.38)
    # Solve for the diameter of small arteriole...
    Diam_sa = fsolve(Tension2,Diam_la,args=params)
    k = Tension2(Diam_sa, *params)
    if(Diam_sa<1e-6 or abs(k)>1e-2):
        n = 1
        while(True): # Loop to calculate an optimum diameter to reduce the fluctuations....
            Diam_sa = (Diameter_sa[-1]*1e-6)+(n*1e-7)
            Diam_sa = fsolve(Tension2,Diam_sa,args=params)
            print('Hitting a hole at a diameter (small) = ', Diam_sa, ' Pressure = ', Pres)
            k = Tension2(Diam_sa, *params)
            print('error is small = ',k)
            if(Diam_sa<3e-5 and abs(k)<0.1):
                break
            n+=1
    # Store the pressures and diameters....
    Pressure_in.append(Pres)
    Diameter_la.append(Diam_la * 1e06)
    Diameter_sa.append(Diam_sa * 1e06)

    Pres+=0.5

k = 0
# Loop to calculate the perfusion...
while(k<len(Diameter_la)):
    gradP_tot = (Pressure_in[k] - 14) * 133.33 # Pressure converted from mmHg to N/m^2 by multiplying with 133.33
    Resistance_la = compartment_resistance(vis_la,l_la,Diameter_la[k]*1e-6,n_la) 
    Resistance_sa = compartment_resistance(vis_sa,l_sa,Diameter_sa[k]*1e-6,n_sa) 
    Resistance_c  = compartment_resistance(vis_c,l_c,Diam_c,n_c) 
    Resistance_lv = compartment_resistance(vis_lv,l_lv,Diam_lv,n_lv) 
    Resistance_sv = compartment_resistance(vis_sv,l_sv,Diam_sv,n_sv) 
    Resistance_a =  compartment_resistance(vis_a,l_a,Diam_a,n_a) 
    Resistance_v =  compartment_resistance(vis_v,l_v,Diam_v,n_v)
    Resistance_total = Resistance_sa + Resistance_la + Resistance_c + Resistance_lv + Resistance_sv
    Resistance_total = Resistance_total + Resistance_a # Calculate the total resistance....
    #print(Resistance_total)
    Q_tot = gradP_tot / Resistance_total
    Diam_a1 = Diam_a + (2 * d_t)
    Diam_c1 = Diam_c + (2 * d_t)
    Diam_la1 = Diam_la + (2 * d_t)
    Diam_sa1 = Diam_sa + (2 * d_t)
    # Calculate the total volume of the arteries and segments...
    vol_a = np.pi * 0.25 * (Diam_a1 * Diam_a1) * l_a * n_a
    vol_la = np.pi * 0.25 * (Diam_la1 * Diam_la1) * l_la * n_la
    vol_sa = np.pi * 0.25 * (Diam_sa1 * Diam_sa1) * l_sa * n_sa
    vol_c = np.pi * 0.25 * (Diam_c1 * Diam_c1) * l_c * n_c
    vol_lv = np.pi * 0.25 * (Diam_lv * Diam_lv) * l_lv * n_lv
    vol_sv = np.pi * 0.25 * (Diam_sv * Diam_sv) * l_sv * n_sv
    vol_v = np.pi * 0.25 * (Diam_v * Diam_v) * l_v * n_v
    vol_tot = vol_la + vol_sa + vol_c + vol_lv + vol_sv + vol_a + vol_v
    perfuse = Q_tot / (vol_tot) # Calculate the total perfusion....
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
for d in csv.DictReader(open('./Carlson_2008(myo+shear+meta).csv')):
    Test_Pressure.append(float(d['Pressure']))
    Test_Diameter.append(float(d['Diameter']))
for d in csv.DictReader(open('./perfusion(myo+shear+meta).csv')):
    Test_Pressure1.append(float(d['Pressure']))
    Test_perfusion.append(float(d['Perfusion']))
# print('Pressure = ', Test_Pressure)
# print('Diameter = ', Test_Diameter)
plt.figure(figsize=(10,9))
plt.subplot(211)
plt.xlabel('Pressure(mmHg)')
plt.ylabel('Diameter(micrometer)')
plt.xlim(0, 200)
plt.ylim(0,180)
plt.plot(Test_Pressure, Test_Diameter, 'b')
plt.plot(Pressure_in,Diameter_la,'r')
plt.subplot(212)
plt.xlabel('Pressure(mmHg)')
plt.ylabel('Perfusion')
plt.xlim(0, 200)
plt.ylim(0,3)
plt.plot(Test_Pressure1, Test_perfusion, 'b')
plt.plot(Pressure_in,perfusion_norm,'r')
plt.show()
