import sys
sys.path.append('./src')
import Perfusion(Myo+Shear+Meta) as perf
def is_float(str):
    try:
        num = float(str)
    except ValueError:
        return False
    return True

def read_file(filename):
    Test_Pressure = []
    Test_Diameter = []
    h = open(filename)
    lines = h.readlines()
    for line in lines:
        k = line.split(',')
        if(is_float(k[0])):
           num = float(k[0])
           num1 = float(k[1])
           Test_Pressure.append(num)
           Test_Diameter.append(num1)
        else:
            continue
    return {'Pressure':Test_Pressure,'Diameter':Test_Diameter}

def check_perfusion_meta():
    Test_Pressure = []
    Test_Diameter = []
    D1= read_file('./files/perfusion(myo+shear+meta).csv')
    for l in D1['Pressure']:
        print(' The values are',l)



