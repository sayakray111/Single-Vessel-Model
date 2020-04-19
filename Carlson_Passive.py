import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
import math
import csv
import scipy.io
import os.path
#from pathlib import Path

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
# The length of the various segments.
l_c = 0.05*0.01
l_v = 0.61*0.01
l_a = 0.61*0.01
l_la = 0.59*0.01
l_lv = 0.59*0.01
l_sv = 0.36*0.01
l_sa = 0.36*0.01
# Numbers of various vessels.
n_c =  5515
n_la = 1
n_a =  1
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
shear  = 0.0
# Function for calculating the tension in the segment... This is fed to fsolve to solve for the diameter at which net Tension is zero.
def Tension2(D,*params):
    P,Cpass,Cpassd,Cact,Cactd,Cactdd,Cmyo,Cshear,Cmeta,Ctoned,N1,D0,vis_1,meta = params
    lam_D = D/D0
    Q = Q_tot/N1
    Q = 0
    shear = (32 * Q * vis_1)/((np.pi) * math.pow(D, 3))
    shear = 5.5 # 55 dyne/cm^2 converted to 5.5 N/m^2
    Stone = (Cmyo*(P*D*0.5))+Ctoned-(Cshear*shear)
    
    Act = 1/(1+np.exp(-Stone))

    p1 = (lam_D-1)*Cpassd

    p2 = np.exp(p1)
    Tpass = Cpass * p2
    s1 = (lam_D-Cactd)/(Cactdd)
    s2 = math.pow(s1,2)
    s3 = np.exp(-s2)
    Tact = Cact * s3
    Ttot = Tpass#+(Act*Tact)
    #print('Activation = ',Act,'Passive Tension = ',Tpass,'Active Tension = ',(Tact*Act))
    #print(Ttot-(P*D*0.5))
    return (Ttot)-(P*D*0.5)

#Viscosity values....
vis_a = 2.22646*0.001 #viscosity of artery (Pa.s), directly from Arciero et al
vis_la = 2.1087*0.001#viscosity of large arteriole (Pa.s), directly from Arciero et al
vis_sa = 3.5542*0.001 #viscosity of small arteriole(Pa.s), directly from Arciero et al
vis_c = 9.0517*0.001 #viscosity of capillary (Pa.s), directly from Arciero et al
vis_sv = 2.5588*0.001 #viscosity of small venule (Pa.s), directly from Arciero et al
vis_lv = 2.3356*0.001  #viscosity of large venule (Pa.s), directly from Arciero et al
vis_v = 2.5007*0.001 #viscosity of vein (Pa.s), directly from Arciero et al

# Control Diameters of the segments.... Non vasoactive segments retain control diameters 
# vasoactive segments change their diameter with pressure..
Diam_la = 0.001 # Initial guess for large arteriole diameter
Diam_sa = 0.001 # Initial guess for small arteriole diameter
Diam_c = 6*1e-6 # 6 micrometer converted to meters Diameter of the capillary
Diam_lv = 119.1*1e-6 # 119.1 micrometers conveted to meters Diameter of the large vein
Diam_sv = 23.5*1e-6 # 23.5 micrometers converted to meters Diameter of the small vein
Diam_lac = 65.2*1e-6 #What is lac? 65.2 micrometers converted to meters.. Control diameter of the large arteriole.
Diam_sac = 14.8*1e-6 #what is sac? 14.8 micrometers converted to meters.. Control diameter of the small arteriole.

#add a comment explaining why this is done
# This is done to calculate the diameter of the artery and the vein from their resistances...
Resistance_a = 1.88*1e13 # Converted to SI units ... 
G = (Resistance_a*np.pi*n_a)/(128*vis_a*l_a)
Diam_a = math.pow((1/G),0.25)
Resistance_v = 0.19*1e13
G = (Resistance_v*np.pi*n_v)/(128*vis_v*l_v)
Diam_v = math.pow((1/G),0.25)

#What is happening here?
# Here the diameters of the thicknesses of the tissue layers around the arteries these values are never used except in the perfusion calculations.
d_t = 18.8 * 1e-6 # Converted to meters from micrometers
Diam_a1 = Diam_a + (2 * d_t)
Diam_lac1 = Diam_lac + (2 * d_t)
Diam_sac1 = Diam_sac + (2 * d_t)
Diam_c1 = Diam_c + (2*d_t)
# Initialise some terms...
perfusion = []
perfuse_100 = 0
Pressure_sa = []
Diameter_sa = []
Diameter_la = []
Pressure_la = []
Pres = 20
D111 = 0.0
perfuse = 0.0
Pressure_in = []
# Coefficients for large arteriole...
Cpass_la = 1042.99*0.001 # Converted from dyne/cm to N/m
Cpassd_la = 8.293 # unitless so remain same
Cact_la = 1596.3*0.001 # Converted from dyne/cm to N/m
Cactd_la = 0.6804 # Unitless
Cactdd_la = 0.2905 # Unitless
Cmyo_la = 0.0101*1000 # Convert cm/dyne to m/N
Cshear_la = 0.0258*10 # Convert cm^2/dyne to m^2/N
Cmeta_la = 30*1e05 # Convert micromole/cm to 1mole/m^3/m
Ctoned_la = -2.22 # unitless so remain same
Ctonedd_la = 10.11 # unitless so remain same
D0_la = 156.49 * 1e-6 # convert micrometer to m
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
    gradP_tot = (Pres-14) * 133.33 # Pressure converted from mmHg to N/m^2 by multiplying with 133.33
    #Calculate resistances via a function, not one by one
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
    Q_tot = gradP_tot/Resistance_total

    ## These are calculations with respect to the large arteriole.... 
    ## The parameters used here to calculate the large arteriole diameter corresponding to a particular large arteriole midpoint pressure
    P1 = (Pres*133.333)-(Q_tot*Resistance_a)
    P2 = (Pres*133.333)-(Q_tot*Resistance_a)-(Q_tot*Resistance_la)
    #P2 = P2 + (Q_tot*Resistance_v)
    ## The mid point pressure for the large arteriole.... (P_la)
    P_la = (P1+P2)*0.5
    D111 = P_la
    
    #define these parameters outside the loop and name them for their artery, make them inputs to your tension function

    ## All the following parameters are for the large arteriole....
    params = (D111,Cpass_la,Cpassd_la,Cact_la,Cactd_la,Cactdd_la,Cmyo_la,Cshear_la,Cmeta_la,Ctoned_la,n_la,D0_la,vis_la,meta_la)
    Diam_la = fsolve(Tension2,Diam_la,args=params)
    k = Tension2(Diam_la,*params)
    # Check if solution has converged if not repeat the problem with that value...
    if(abs(k)>0.00008):
        Diam_la = Diam_la+0.000001
        print(Pres)
        continue
        
    #define these parameters outside the loop and name them for their artery, make them inputs to your tension function
    ## These are calculations with respect to the small arteriole.... 
    ## The parameters used here to calculate the small arteriole diameter corresponding to a particular small arteriole midpoint pressure
    #parameters for the small arteriole...
    # Here the midpoint pressure is being calculated....
    P22 = (Pres*133.333)-(Q_tot*Resistance_a)-(Q_tot*Resistance_la)-(Q_tot*Resistance_sa) # all pressures are in N/m^2
    P11 = (Pres*133.333)-(Q_tot*Resistance_la)-(Q_tot*Resistance_a)
    P_sa = (P11+P22)*0.5
    
    D111 = P_sa
    ## The mid point pressure for the small arteriole.... (P_sa)
    # These parameters are for the small arteriole...
    params = (D111,Cpass_sa,Cpassd_sa,Cact_sa,Cactd_sa,Cactdd_sa,Cmyo_sa,Cshear_sa,Cmeta_sa,Ctoned_sa,n_sa,D0_sa,vis_sa,meta_sa)
    ## Solve for the small arteriole diameter....
    Diam_sa = fsolve(Tension2,Diam_sa,args=params)
    
    #what is this loop supposed to do?
    k = Tension2(Diam_sa,*params)
    # Check if solution has converged if not repeat the problem with that value...
    if(abs(k)>0.00008):
        Diam_sa = Diam_sa+0.000001
        print(Pres)
        continue
        
    ## Store the values of the input pressure and the Diameters of the large and small arterioles.
    Pressure_in.append(Pres)
    Pressure_la.append(P_la/133.33)
    Diameter_la.append(Diam_la * 1e06)
    Pressure_sa.append(P_sa/133.33)
    Diameter_sa.append(Diam_sa*1e06)
    
    

    #why redefine d_t every loop? move calculation of all updated diameters here
    #I do not recall you discussing this at all in your description of this model

    # Here the diameters have been converted to their total values which considers the amount of tissue covering the vessel..
    Diam_la1 = Diam_la +(2*d_t)
    Diam_sa1 = Diam_sa+(2*d_t)
    ## Here the perfusion is calculated on the basis of the formulae given in the report and storing the values...
    vol_a = np.pi * 0.25 * (Diam_a1 * Diam_a1) * l_a * n_a
    vol_la = np.pi * 0.25 * (Diam_la1 * Diam_la1) * l_la * n_la
    vol_sa = np.pi * 0.25 * (Diam_sa1 * Diam_sa1) * l_sa * n_sa
    vol_c = np.pi * 0.25 * (Diam_c1 * Diam_c1) * l_c * n_c
    vol_lv = np.pi * 0.25 * (Diam_lv * Diam_lv) * l_lv * n_lv
    vol_sv = np.pi * 0.25 * (Diam_sv * Diam_sv) * l_sv * n_sv
    vol_v = np.pi * 0.25 * (Diam_v * Diam_v) * l_v * n_v
    vol_tot = vol_la + vol_sa + vol_c + vol_lv + vol_sv + vol_a + vol_v
    perfuse = Q_tot / (vol_tot)
    perfusion.append(perfuse) # Here perfusion is in 1/sec. but the value in the paper is in 1/100 * 60 secs so should be divided by 6000
    
    
    # Calculating the normalised perfusion
    if(Pres==47):
        perfuse_100 = perfuse
        vol_100 = vol_tot
        Resistance_total100 = Resistance_total
        
    Pres+=1
    
    

#Normalising the perfusion values...
perfusion_norm = [(k/perfuse_100) for k in perfusion]
Test_Pressure = []
Test_Pressure1=[]
Test_Diameter = []
Test_perfusion= []
gradP_100 = ((70.9/6000)*(Resistance_total100)*vol_100)/133
print(gradP_100)
for d in csv.DictReader(open('./perfusion(passive).csv')):
    Test_Pressure.append(float(d['Pressure']))
    Test_perfusion.append(float(d['Perfusion']))
for d1 in csv.DictReader(open('./Carlson(2008)_dataP.csv')):
    Test_Pressure1.append(float(d1['Pressure']))
    Test_Diameter.append(float(d1['Diameter']))

plt.figure(figsize=(10,9))
plt.subplot(211)
plt.xlim(0,200)
plt.ylim(0,3)
plt.xlabel('Pressure(mmHg)')
plt.ylabel('perfusion')
plt.plot(Test_Pressure,Test_perfusion,'b')
plt.plot(Pressure_in,perfusion_norm,'r')
plt.subplot(212)
plt.xlim(0,200)
plt.ylim(0,180)
plt.xlabel('Pressure(mmHg)')
plt.ylabel('Diameter')
plt.plot(Test_Pressure1,Test_Diameter,'b')
plt.plot(Pressure_la,Diameter_la,'r')
plt.show()
