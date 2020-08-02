import os
import test_SVM as cr
if __name__ == '__main__':

    export_directory = 'output'
    if not os.path.exists(export_directory):
        os.makedirs(export_directory)

    everymodule_worked = True
    current_dir = os.getcwd()
    Flag1 = cr.check_perfusion1
    Flag2 = cr.check_perfusion2
    Flag3 = cr.check_perfusion3
    Flag4 = cr.check_perfusion4


    if(Flag1 and Flag2 and Flag3 and Flag4):
        everymodule_worked = True
    else:
        everymodule_worked = False
    if(everymodule_worked==True):
        print('The codes have all worked.')
    elif(Flag1==False):
        print('The passive code hasnt worked')
    elif(Flag2==False):
        print('The Myogenic code hasnt worked')
    elif(Flag3==False):
        print('The shear code hasnt worked')
    else:
        print('The Meta code hasnt worked')
