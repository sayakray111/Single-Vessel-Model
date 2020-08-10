import src.Perfusion1 as p1
import src.Perfusion2 as p2
import src.Perfusion3 as p3
import src.Perfusion4 as p4
import src.Perfusion5 as p5
import matplotlib.pyplot as plt
import numpy as np
import csv
def is_float(str):
    try:
        num = float(str)
    except ValueError:
        return False
    return True

def read_file(filename):
    Test_Pressure = []
    Test_Perfusion = []
    h = open(filename)
    lines = h.readlines()
    for line in lines:
        k = line.split(',')
        if(is_float(k[0])):
           num = float(k[0])
           num1 = float(k[1])
           Test_Pressure.append(num)
           Test_Perfusion.append(num1)
        else:
            continue
    return {'Pressure':Test_Pressure,'Perfusion':Test_Perfusion}

def check_perfusion1():
    Test_Pressure = []
    Test_Perfusion = []
    D1= read_file('./files/perfusion(passive).csv')
    for l,k in zip(D1['Pressure'],D1['Perfusion']):
        Test_Pressure.append(l)
        Test_Perfusion.append(k)
    A1 = p1.Perfusion_pass()
    D2 = A1.Calculate()
    Pressure_in = D2['Pressure']
    perfusion_norm = D2['Normalised Perfusion']
    perfusion_norm1 = np.interp(Pressure_in,Test_Pressure,Test_Perfusion).tolist()
    error = 0.0
    k1 = 0
    while(k1<len(perfusion_norm1)):
        error = error+(abs(perfusion_norm1[k1]-perfusion_norm[k1])/perfusion_norm1[k1])
        if(perfusion_norm[k1]>3.0):
            break
        k1+=1
    avg = error/len(perfusion_norm1)
    print('The approximate error is around = ',avg*100)
    if((avg*100)<=1.5):
        return True
    else:
        return False

def check_perfusion2():
    Test_Pressure = []
    Test_Perfusion = []
    D1= read_file('./files/perfusion(myo).csv')
    for l,k in zip(D1['Pressure'],D1['Perfusion']):
        Test_Pressure.append(l)
        Test_Perfusion.append(k)
    D2 = p2.Perfusion_myo()
    Pressure_in = D2['Pressure']
    perfusion_norm = D2['Normalised Perfusion']
    perfusion_norm1 = np.interp(Pressure_in,Test_Pressure,Test_Perfusion).tolist()
    error = 0.0
    k1 = 0
    while (k1 < len(perfusion_norm1)):
        error = error + (abs(perfusion_norm1[k1] - perfusion_norm[k1]) / perfusion_norm1[k1])
        if (perfusion_norm[k1] > max(Test_Perfusion)):
            break
        k1 += 1
    avg = error / len(perfusion_norm1)
    print('The approximate error is around = ', avg * 100)
    if ((avg * 100) <= 1.5):
        return True
    else:
        return False

def check_perfusion3():
    Test_Pressure = []
    Test_Perfusion = []
    D1= read_file('./files/perfusion(myo+shear).csv')
    for l,k in zip(D1['Pressure'],D1['Perfusion']):
        Test_Pressure.append(l)
        Test_Perfusion.append(k)
    D2 = p3.Perfusion_shear()
    Pressure_in = []
    perfusion_norm = []
    Pressure_in = D2['Pressure']
    perfusion_norm = D2['Normalised Perfusion']
    perfusion_norm1 = np.interp(Pressure_in,Test_Pressure,Test_Perfusion).tolist()
    error = 0.0
    k1 = 0
    while (k1 < len(perfusion_norm1)):
        error = error + (abs(perfusion_norm1[k1] - perfusion_norm[k1]) / perfusion_norm1[k1])
        if (perfusion_norm[k1] > max(Test_Perfusion)):
            break
        k1 += 1
    avg = error / len(perfusion_norm1)
    print('The approximate error is around = ', avg * 100)
    if ((avg * 100) <= 1.5):
        return True
    else:
        return False

def check_perfusion4():
    Test_Pressure = []
    Test_Perfusion = []
    D1= read_file('./files/perfusion(myo+shear+meta).csv')
    for l,k in zip(D1['Pressure'],D1['Perfusion']):
        Test_Pressure.append(l)
        Test_Perfusion.append(k)
    D2 = p4.Perfusion_meta()
    Pressure_in = []
    perfusion_norm = []
    Pressure_in = D2['Pressure']
    perfusion_norm = D2['Normalised Perfusion']
    perfusion_norm1 = np.interp(Pressure_in,Test_Pressure,Test_Perfusion).tolist()
    error = 0.0
    k1 = 0
    while (k1 < len(perfusion_norm1)):
        error = error + (abs(perfusion_norm1[k1] - perfusion_norm[k1]) / perfusion_norm1[k1])
        if (perfusion_norm[k1] > max(Test_Perfusion)):
            break
        k1 += 1
    avg = error / len(perfusion_norm1)
    print('The approximate error is around = ', avg * 100)
    if ((avg * 100) <= 1.5):
        return True
    else:
        return False

def plot_perfusion1():
    Test_Pressure = []
    Test_Pressure1 = []
    Test_Pressure2 = []
    Test_Diameter = []
    Test_Diameter1 = []
    Test_perfusion = []
    # gradP_100 = ((70.9 / 6000) * (Resistance_total100) * vol_100) / 133
    # print(Diameter_sa)
    # print(S3)
    F1 = p1.Perfusion_pass()
    F = F1.Calculate()
    for d in csv.DictReader(open('./files/Carlson(2008)_dataP.csv')):
        Test_Pressure.append(float(d['Pressure']))
        Test_Diameter.append(float(d['Diameter']))
    for d in csv.DictReader(open('./files/perfusion(passive).csv')):
        Test_Pressure1.append(float(d['Pressure']))
        Test_perfusion.append(float(d['Perfusion']))
    Pressure_la = F['Pressure(LA)']
    Diameter_la = F['Diameter(LA)']
    Pressure_in = F['Pressure']
    perfusion = F['Perfusion']
    # print('Pressure = ', Test_Pressure)
    # print('Diameter = ', Test_Diameter)
    plt.figure(figsize=(10, 9))
    plt.subplot(311)
    plt.xlabel('Pressure(mmHg)')
    plt.ylabel('Diameter(micrometer)')
    plt.xlim(0, 200)
    plt.ylim(0, 180)
    plt.plot(Test_Pressure, Test_Diameter, 'b')
    plt.plot(Pressure_la, Diameter_la, 'r')
    plt.subplot(312)
    plt.xlabel('Pressure(mmHg)')
    plt.ylabel('Perfusion')
    plt.xlim(0, 200)
    plt.ylim(0, 3)
    plt.plot(Test_Pressure1, Test_perfusion, 'b')
    plt.plot(Pressure_in, perfusion, 'r')
    plt.show()







