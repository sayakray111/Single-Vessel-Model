import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
import math
import csv
import scipy.io
import os.path
from pathlib import Path
def compartment_resistance(viscosity,length,diameter,number_in_generation):
    resistance  = (128. * viscosity* length) / (np.pi * math.pow(diameter, 4.) * number_in_generation)
    return resistance
mat = scipy.io.loadmat('./DAnometa2vals.mat')
Pvc = float(mat.get('Pvc') / 1333)
Pac = float(mat.get('Pac') / 1333)
#calculate the value of Q to be used for the evaluation of shear stress
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
Q_sa = Q_tot/143
Q_lv = Q_tot
Q_sv = Q_tot/143
Q_c = Q_tot/5515


# The length of the various segments.
l_c = 0.05*0.01
l_v = 0.61*0.01
l_a = 0.61*0.01
l_la = 0.59*0.01
l_lv = 0.59*0.01
l_sv = 0.36*0.01
l_sa = 0.36*0.01
# the pressure drop in the segments.
P_la = 0
P_lv = 0
P_sa = 0
P_sv = 0
P_c = 0
# Numbers of various vessels.
n_c =  5515
n_la = 1
n_a =  1
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
shear  = 0

#P = 100*133.333
Diameter = []
Pressure = []
f1  = []
f2  = []
Tes = []
vis = 7*1e-4

sign = lambda x: math.copysign(1, x)
N1 = 0

def Tension2(D,*params):
    
    P,Cpass,Cpassd,Cact,Cactd,Cactdd,Cmyo,Cshear,Cmeta,Ctoned,N1,D0,vis_1,meta = params
    lam_D = D/D0
    Q = Q_tot/N1
    Q = 0
    shear = (32 * Q * vis_1)/((np.pi) * math.pow(D, 3))
    shear = 5.5
    Stone = (Cmyo*(P*D*0.5))+Ctoned-(Cshear*shear)
    
    Act = 1/(1+np.exp(-Stone))

    p1 = (lam_D-1)*Cpassd

    p2 = np.exp(p1)
    Tpass = Cpass * p2
    s1 = (lam_D-Cactd)/(Cactdd)
    s2 = math.pow(s1,2)
    s3 = np.exp(-s2)
    Tact = Cact * s3
    Ttot = Tpass+(Act*Tact)
    #print('Activation = ',Act,'Passive Tension = ',Tpass,'Active Tension = ',(Tact*Act))
    #print(Ttot-(P*D*0.5))
    return (Ttot)-(P*D*0.5)



def Activation(T):
    Stim = (Cmyo*T)-Ctoned
    return 1/(1+np.exp(-Stim))
#Viscosity values....
vis_a = 0.022646*0.1
vis_la = 0.021087*0.1
vis_sa = 0.035542*0.1
vis_c = 0.090517*0.1
vis_sv = 0.025588*0.1
vis_lv = 0.023356*0.1
vis_v = 0.025007*0.1
# Control Diameters of the segments....
Diam_la = 0.001
Diam_sa = 0.001
Diam_c = 6*1e-6
Diam_lv = 119.1*1e-6
Diam_sv = 23.5*1e-6
Diam_lac = 65.2*1e-6
Diam_sac = 14.8*1e-6
d_t = 18.8 * 1e-6
Resistance_a = 1.88*1e13
G = (Resistance_a*np.pi*n_a)/(128*vis_a*l_a)
Diam_a = math.pow((1/G),0.25)
Resistance_v = 0.19*1e13
G = (Resistance_v*np.pi*n_v)/(128*vis_v*l_v)
Diam_v = math.pow((1/G),0.25)
Diam_a1 = Diam_a + (2 * d_t)
Diam_lac1 = Diam_lac + (2 * d_t)
Diam_sac1 = Diam_sac + (2 * d_t)
perfuse = 1
perfusion = []
perfuse_100 = 0
Pressure_sa = []
Diameter_sa = []
Diameter_la = []
Pressure_la = []
Pres = 14
D111 = 0.0
perfuse = 0.0
Pres = 20
Pressure_in = []
# Coeffiecients for large arteriole...
Cpass_la = 1.043
Cpassd_la = 8.293
Cact_la = 1.596
Cactd_la = 0.6804
Cactdd_la = 0.2905
Cmyo_la = 10.1
Cshear_la = 0.258
Cmeta_la = 3*1e6
Ctoned_la = -2.22
Ctonedd_la = 10.11
D0_la = 156.49 * 1e-6
# Coefficients for the small arteriole...
Cpass_sa = 0.2599
Cpassd_sa = 11.467
Cact_sa = 0.274193
Cactd_sa = 0.750
Cactdd_sa = 0.383544
Cmyo_sa = 35.867153
Cshear_sa = 0.258
Cmeta_sa = 3*1e06
Ctoned_sa = -0.53
Ctonedd_sa = 10.66
D0_sa = 38.99*1e-6
# meta calculated in the control state....
Tlac = 100*133.33*Diam_lac*0.5
meta_la = ((Cmyo_la*Tlac)-(Cshear_la*5.5)+Ctonedd_la)/Cmeta_la
print(meta_la)
Tsac = 100*133.33*Diam_sac*0.5
meta_sa = ((Cmyo_sa*Tsac)-(Cshear_sa*5.5)+Ctonedd_sa)/Cmeta_sa
print(meta_sa)
while(Pres<=200):
    #large arteriole...
    gradP_tot = (Pres-14) * 133.33
    Resistance_la  = compartment_resistance(vis_la,l_la,Diam_la,n_la) 
    Resistance_sa  = compartment_resistance(vis_sa,l_sa,Diam_sa,n_sa) 
    Resistance_c   = compartment_resistance(vis_c,l_c,Diam_c,n_c) 
    Resistance_lv  = compartment_resistance(vis_lv,l_lv,Diam_lv,n_lv) 
    Resistance_sv  = compartment_resistance(vis_sv,l_sv,Diam_sv,n_sv) 
    Resistance_lac = compartment_resistance(vis_la,l_la,Diam_lac,n_la) 
    Resistance_sac = compartment_resistance(vis_sa,l_sa,Diam_sac,n_sa) 
    Resistance_a = compartment_resistance(vis_a,l_a,Diam_a,n_a) 
    Resistance_v = compartment_resistance(vis_v,l_v,Diam_v,n_v) 
    Resistance_total = Resistance_sa + Resistance_la + Resistance_c + Resistance_lv + Resistance_sv
    Resistance_totalc = Resistance_sac + Resistance_lac + Resistance_c + Resistance_lv + Resistance_sv
    Resistance_total = Resistance_total + Resistance_a + Resistance_v
    #Resistance_totalc = Resistance_totalc + Resistance_a + Resistance_v
    Q_tot = gradP_tot/Resistance_total
    #Rat = ((Pres-Pvc)/(Pac-Pvc))
    #
    P1 = (Pres*133.333)-(Q_tot*Resistance_a)
    P2 = (Pres*133.333)-(Q_tot*Resistance_a)-(Q_tot*Resistance_la)
    #P2 = P2 + (Q_tot*Resistance_v)
    P_la = (P1+P2)*0.5
    D111 = P_la
    
    #X0 = Dekkers(Tension2,0,0.00005,D111)
    params = (D111,Cpass_la,Cpassd_la,Cact_la,Cactd_la,Cactdd_la,Cmyo_la,Cshear_la,Cmeta_la,Ctoned_la,n_la,D0_la,vis_la,meta_la)
    Diam_la = fsolve(Tension2,Diam_la,args=params)
    
    k = Tension2(Diam_la,*params)
    #print(k)
    if(abs(k)>0.00008):
        Diam_la = Diam_la+0.000001
        print(Pres)
        continue
    #small arteriole...
    P22 = (Pres*133.333)-(Q_tot*Resistance_a)-(Q_tot*Resistance_la)-(Q_tot*Resistance_sa)
    P11 = (Pres*133.333)-(Q_tot*Resistance_la)-(Q_tot*Resistance_a)
    P_sa = (P11+P22)*0.5
    D111 = P_sa
    #X0 = Dekkers(Tension2,0,0.00005,D111)
    params = (D111,Cpass_sa,Cpassd_sa,Cact_sa,Cactd_sa,Cactdd_sa,Cmyo_sa,Cshear_sa,Cmeta_sa,Ctoned_sa,n_sa,D0_sa,vis_sa,meta_sa)
    Diam_sa = fsolve(Tension2,Diam_sa,args=params)
    k = Tension2(Diam_sa,*params)
    if(abs(k)>0.00008):
        Diam_sa = Diam_sa+0.000001
        print(Pres)
        continue
    Pressure_in.append(Pres)
    Pressure_la.append(P_la/133.33)
    Diameter_la.append(Diam_la * 1e06)
    Pressure_sa.append(P_sa/133.33)
    Diameter_sa.append(Diam_sa*1e06)
    d_t = 18.8*1e-6
    
    Diam_a1 = Diam_a+(2*d_t)
    Diam_la1 = Diam_la+(2*d_t)
    Diam_sa1 = Diam_sa+(2*d_t)
    Diam_lac1 = Diam_lac + (2*d_t)
    Diam_sac1 = Diam_sac + (2*d_t)
    Diam_c1 = Diam_c + (2*d_t)
    vol_a = np.pi * 0.25 * (Diam_a1 * Diam_a1) * l_a * n_a
    vol_la = np.pi * 0.25 * (Diam_la1 * Diam_la1) * l_la * n_la
    vol_sa = np.pi * 0.25 * (Diam_sa1 * Diam_sa1) * l_sa * n_sa
    vol_c = np.pi * 0.25 * (Diam_c1 * Diam_c1) * l_c * n_c
    vol_lv = np.pi * 0.25 * (Diam_lv * Diam_lv) * l_lv * n_lv
    vol_sv = np.pi * 0.25 * (Diam_sv * Diam_sv) * l_sv * n_sv
    vol_v = np.pi * 0.25 * (Diam_v * Diam_v) * l_v * n_v
    vol_tot = vol_la + vol_sa + vol_c + vol_lv + vol_sv + vol_a + vol_v
    perfuse = Q_tot / (vol_tot)
    perfusion.append(perfuse/6000)
    if(Pres==100):
        perfuse_100 = perfuse/6000
        vol_100 = vol_tot
        Resistance_total100 = Resistance_total
    Pres+=1
    
    

perfusion_norm = [(k/perfuse_100) for k in perfusion]
Test_Pressure = []
Test_Diameter = []
Test_perfusion= []
gradP_100 = ((70.9/6000)*(Resistance_total100)*vol_100)/133
print(gradP_100)
for d in csv.DictReader(open('./perfusion(myo).csv')):
    Test_Pressure.append(float(d['Pressure']))
    Test_perfusion.append(float(d['Perfusion']))

print('Pressure = ', Test_Pressure)
print('Diameter = ', Test_Diameter)
plt.xlim(0,200)
plt.ylim(0,3)
plt.xlabel('Pressure(mmHg)')
plt.ylabel('perfusion(normalised)')
plt.plot(Test_Pressure,Test_perfusion,'b')
plt.plot(Pressure_in,perfusion_norm,'r')
plt.show()
