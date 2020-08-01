
import src.Perfusion1 as p1
import src.Perfusion2 as p2
import src.Perfusion3 as p3
import src.Perfusion4 as p4
import src.Perfusion5 as p5
import matplotlib.pyplot as plt
import numpy as np

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
    error = sum(abs(perfusion_norm[i]-perfusion_norm1[i])/perfusion_norm1[i] for i in range(0,len(perfusion_norm1)))
    error = error/len(perfusion_norm1)
    print('The approximate error is around = ',error*100)


check_perfusion4()





