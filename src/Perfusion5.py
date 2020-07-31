# This code derives the Pressure diameter curves in the figure for three different flow rates. 
import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
import math
import csv
def Perfusion5():
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
    vis = 7e-4
    Q = (3 / 3) * 1e-09
    sign = lambda x: math.copysign(1, x)

    def Tension2(D, P):
        lam_D = D / D0
        shear = (32 * Q * vis) / ((np.pi) * math.pow(D, 3))
        Stone = (Cmyo * (P * D * 0.5)) + Ctoned - (Cshear * shear)
        Act = 1 / (1 + np.exp(-Stone))
        p1 = (lam_D - 1) * Cpassd
        # print(p1)
        p2 = np.exp(p1)
        Tpass = Cpass * p2
        s1 = (lam_D - Cactd) / (Cactdd)
        s2 = math.pow(s1, 2)
        s3 = np.exp(-s2)
        Tact = Cact * s3
        Ttot = Tpass + (Act * Tact)
        # print('Activation = ',Act,'Passive Tension = ',Tpass,'Active Tension = ',(Tact*Act))
        # print(Ttot-(P*D*0.5))
        return (Ttot) - (P * D * 0.5)

    def Activation(T):
        Stim = (Cmyo * T) - Ctoned
        return 1 / (1 + np.exp(-Stim))

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
        # print(k)
        if (abs(k) > 0.0008):
            X0 = Diam + 0.000001

            continue

        Pressure.append(Pres)
        Diameter.append(Diam * 1e06)
        Pres += 1

    Test_Pressure = []
    Test_Diameter = []

    for d in csv.DictReader(open('D:/Photos/My_Work/Carlson(2008)_data60l.csv')):
        Test_Pressure.append(float(d['Pressure']))
        Test_Diameter.append(float(d['Diameter']))

    print('Pressure = ', Test_Pressure)
    print('Diameter = ', Test_Diameter)
    Stimuli = [Activation(ki) for ki in Tes]
    plt.xlim(0, 120)
    plt.ylim(0, 180)

    plt.xlabel('Pressure')
    plt.ylabel('Diameter')
    plt.plot(Test_Pressure, Test_Diameter)
    plt.plot(Pressure, Diameter)
    plt.show()



        
        