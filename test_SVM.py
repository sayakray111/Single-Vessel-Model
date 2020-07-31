
import src.Perfusion1
import src.Perfusion2
import src.Perfusion3
import src.Perfusion4
import src.Perfusion5
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

def check_perfusion4():
    Test_Pressure = []
    Test_Perfusion = []
    D1= read_file('../files/perfusion(myo+shear+meta).csv')
    for l,k in zip(D1['Pressure'],D1['Perfusion']):
        Test_Pressure.append(l)
        Test_Perfusion.append(k)






