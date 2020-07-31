import sys
sys.path.append('/hpc/sray036/Single-Vessel-Model/src')

import matplotlib.pyplot as plt

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

def check_perfusion_meta():
    Test_Pressure = []
    Test_Perfusion = []
    D1= read_file('../files/perfusion(myo+shear+meta).csv')
    for l,k in zip(D1['Pressure'],D1['Perfusion']):
        Test_Pressure.append(l)
        Test_Perfusion.append(k)
    plt.xlim(0, 200)
    plt.ylim(0, 3)
    plt.xlabel('Pressure')
    plt.ylabel('Perfusion')
    plt.plot(Test_Pressure, Test_Perfusion,'r')
    plt.show()





