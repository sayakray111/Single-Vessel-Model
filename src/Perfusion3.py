# Written by Sayak Ray of the university of Auckland
# This code creates the perfusion-pressure curve for myo
import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
import math
import csv
import scipy.io
import scipy.integrate as integrate

class Perfusion_shear():
    def __init__(self):
        self.perfusion_norm = []
        self.Diameter_sa = []
        self.perfusion = []
        self.Diameter_la = []
        self.Pressure_in = []
        self.Pressure_la = []
        self.Pressure_sa = []
        #Length of the segments....
        self.l_c = 0.05 * 0.01
        self.l_v = 0.61 * 0.01
        self.l_a = 0.61 * 0.01
        self.l_la = 0.59 * 0.01
        self.l_lv = 0.59 * 0.01
        self.l_sv = 0.36 * 0.01
        self.l_sa = 0.36 * 0.01
        # Numbers of various vessels adapted from the paper....
        self.n_c = 5515
        self.n_la = 1
        self.n_a = 1
        self.n_v = 1
        self.n_sa = 143
        self.n_lv = 1
        self.n_sv = 143
        # volume of the segments...
        self.vol_la = 0.0
        self.vol_sa = 0.0
        self.vol_c = 0.0
        self.vol_sv = 0.0
        self.vol_lv = 0.0
        self.vol_tot = 0.0
        self.shear = 0.0
        self.meta = 0.0
        # Viscosities of the segments...
        self.vis_a = 0.0226 * 0.1  # viscosity of artery (Pa.s), directly from Arciero et al
        self.vis_la = 0.0211 * 0.1  # viscosity of large arteriole (Pa.s), directly from Arciero et al
        self.vis_sa = 0.0355 * 0.1  # viscosity of small arteriole(Pa.s), directly from Arciero et al
        self.vis_c = 0.0905 * 0.1  # viscosity of capillary (Pa.s), directly from Arciero et al
        self.vis_sv = 0.0256 * 0.1  # viscosity of small venule (Pa.s), directly from Arciero et al
        self.vis_lv = 0.0234 * 0.1  # viscosity of large venule (Pa.s), directly from Arciero et al
        self.vis_v = 0.0250 * 0.1  # viscosity of vein (Pa.s), directly from Arciero et al
        # Coefficients for metabolic response... small arteriole. large arteriole... adapted from the paper.
        self.d_tissue = 18.8 * 1e-6
        self.M_0 = (8.28 / 60) * 0.01
        self.c0 = 0.5
        self.H_D = 0.4
        self.H_T = 0.3
        self.R_0 = 1.4e-3
        self.R_1 = 0.891
        self.zeta_H = 2.7
        self.P_50 = 26.0 * 133.322
        self.S_0 = 0.97
        self.C_0_la = 0.5 * 1e-3
        self.k_d = 2.0 * 1e-6
        self.tau_d = 1
        self.tau_a = 60
        self.L_0 = 1.0e-2
        self.x0_la = 0.61e-2
        self.x0_sa = 1.2e-2
        self.x0_c = 1.56e-2
        self.x0_sv = 1.61e-2
        self.x0_lv = 1.97e-2
        self.x0_v = 2.56e-2
        self.xmp_la = 0.90e-2
        self.xmp_sa = 1.38e-2
        self.xend_la = 1.2e-2
        self.xend_sa = 1.56e-2
        self.xend_c = 1.61e-2
        self.xend_sv = 1.97e-2
        self.xend_lv = 2.56e-2
        self.xD_la = 0.90
        self.xD_sa = 1.38
        self.d_t = 18.8 * 1e-6
        # Coefficients for large arteriole.... adapted from the paper.
        self.Cpass_la = 1.043
        self.Cpassd_la = 8.293
        self.Cact_la = 1.596
        self.Cactd_la = 0.6804
        self.Cactdd_la = 0.2905
        self.Cmyo_la = 10.1
        self.Cshear_la = 0.258
        self.Cmeta_la = 3 * 1e06
        self.Ctoned_la = -5.22
        self.Ctonedd_la = 10.11
        self.D0_la = 156.49 * 1e-6
        # Coefficients for small arteriole... adapted from the paper.
        self.Cpass_sa = 0.2599
        self.Cpassd_sa = 11.467
        self.Cact_sa = 0.274193
        self.Cactd_sa = 0.750
        self.Cactdd_sa = 0.384
        self.Cmyo_sa = 35.9
        self.Cshear_sa = 0.258
        self.Ctoned_sa = -15.3
        self.Cmeta_sa = 3 * 1e06
        self.D0_sa = 38.99 * 1e-6
        self.Ctonedd_sa = 10.66
        # This is the coefficients for the arterioles
        # Diameters of the artery.....This is done to calculate the diameter of the artery and the vein from their resistances...
        Resistance_a = 1.88 * 1e11 * 133
        G = (128 * self.vis_a * self.l_a) / (np.pi * Resistance_a * self.n_a)
        self.Diam_a = math.pow(G, 0.25)
        print('Diameter of the main artery is = ', self.Diam_a * 1e06)
        # Diameters of the vein....This is done to calculate the diameter of the artery and the vein from their resistances...
        Resistance_v = 0.19 * 1e11 * 133
        G = (128 * self.vis_v * self.l_v) / (np.pi * Resistance_v * self.n_v)
        self.Diam_v = math.pow(G, 0.25)
        print('Diameter of the main vein is = ', self.Diam_v * 1e06)
        self.Diam_c = 6 * 1e-6  # Diameter of the capillary adapted from the paper
        self.Diam_lv = 119.1 * 1e-6  # Diameter of the large venule adapted from the paper
        self.Diam_sv = 23.5 * 1e-6  # Diameter of the small venule adapted from the paper

        self.d_t = 18.8 * 1e-6  # layer of tissue around the vessels...
        # Control Diameters of the large arteriole adapted from the paper
        self.Diam_lac = 65.2 * 1e-6
        # Control Diameters of the small arteriole adapted from the paper.
        self.Diam_sac = 14.8 * 1e-6
        # Here the diameters of the thicknesses of the tissue layers around the arteries are added to the original diameters..
        # these values are never used except in the perfusion calculations.
        #self.Diam_la1 = self.Diam_la + (2 * self.d_t)
        self.Diam_lac1 = self.Diam_lac + (2 * self.d_t)
        #self.Diam_sa1 = self.Diam_sa + (2 * self.d_t)
        self.Diam_sac1 = self.Diam_sac + (2 * self.d_t)
        self.Diam_c1 = self.Diam_c + (2 * self.d_t)
        self.Diam_sv1 = self.Diam_sv + (2 * self.d_t)
        self.Diam_lv1 = self.Diam_lv + (2 * self.d_t)

    def compartment_resistance(self,viscosity, length, diameter, number_in_generation):
        resistance = (128. * viscosity * length) / (np.pi * math.pow(diameter, 4.) * number_in_generation)
        return resistance



    # This function is input into fsolve to calculate the Optimum Diameter for a particular pressure....
    def Tension2(self,Diam,Pres):
        Diam_la = Diam[0]
        Diam_sa = Diam[1]
        # Determine which diameter is which....
        # Calculate the Discharge for the given diameters and their initial guesses...
        gradP_tot = (Pres - 12.91) * 133.33  # Pressure converted from mmHg to N/m^2 by multiplying with 133.33
        # End vein pressure assumed to be 14mmHg
        # Calculating the resistances here ....
        Resistance_la = self.compartment_resistance(self.vis_la, self.l_la, Diam_la, self.n_la)
        Resistance_sa = self.compartment_resistance(self.vis_sa, self.l_sa, Diam_sa, self.n_sa)
        Resistance_c = self.compartment_resistance(self.vis_c, self.l_c, self.Diam_c, self.n_c)
        Resistance_lv = self.compartment_resistance(self.vis_lv, self.l_lv, self.Diam_lv, self.n_lv)
        Resistance_sv = self.compartment_resistance(self.vis_sv, self.l_sv, self.Diam_sv, self.n_sv)
        Resistance_lac = self.compartment_resistance(self.vis_la, self.l_la, self.Diam_lac, self.n_la)
        Resistance_sac = self.compartment_resistance(self.vis_sa, self.l_sa, self.Diam_sac, self.n_sa)
        Resistance_a = self.compartment_resistance(self.vis_a, self.l_a, self.Diam_a, self.n_a)
        Resistance_v = self.compartment_resistance(self.vis_v, self.l_v, self.Diam_v, self.n_v)
        Resistance_total = Resistance_sa + Resistance_la + Resistance_c + Resistance_lv + Resistance_sv
        Resistance_totalc = Resistance_sac + Resistance_lac + Resistance_c + Resistance_lv + Resistance_sv
        Resistance_total = Resistance_total + Resistance_a  # The total resistance is being calculated...
        Q_tot = gradP_tot / Resistance_total  # Calculates the value of total discharge...

        # -----------------------------------------------------
        Q_la = Q_tot / self.n_la  # Calculate the discharge for the specific case...
        Q_sa = Q_tot / self.n_sa
        shear_la = (32 * Q_la * self.vis_la) / ((np.pi) * math.pow(Diam_la, 3))
        P1 = (Pres * 133.333) - (Q_tot * Resistance_a)
        P2 = (Pres * 133.333) - (Q_tot * Resistance_a) - (Q_tot * Resistance_la)
        Pmid_la = (P1 + P2) * 0.5  # The midpoint diameter of the large arteriole...

        shear_sa = (32 * Q_sa * self.vis_sa) / ((np.pi) * math.pow(Diam_sa, 3))
        P22 = (Pres * 133.333) - (Q_tot * Resistance_a) - (Q_tot * Resistance_la) - (Q_tot * Resistance_sa)
        P11 = (Pres * 133.333) - (Q_tot * Resistance_la) - (Q_tot * Resistance_a)
        Pmid_sa = (P11 + P22) * 0.5  # The midpoint diameter of the small arteriole...
        # shear_la = 5.5
        # shear_sa = 5.5
        # shear_la = 0.0
        # shear_sa = 0.0
        # Calculating the consumption of ATP...
        cons = (Diam_la, Diam_sa, Q_tot)
        # Calculating the metabolic constant (SCR) ...
        # meta_la = SCR(xD_la, *cons)
        # meta_sa = SCR(xD_sa, *cons)
        # meta_la = 4.360341193333333e-06
        # meta_sa = 4.261014926666667e-06
        lam_D_la = Diam_la / self.D0_la  # The ratio of diameters...
        lam_D_sa = Diam_sa / self.D0_sa
        Stone_la = (self.Cmyo_la * (Pmid_la * Diam_la * 0.5)) + self.Ctoned_la - (self.Cshear_la * shear_la)  # - (Cmeta_la * meta_la)  # The stone is the tone of the vessel...
        Stone_sa = (self.Cmyo_sa * (Pmid_sa * Diam_sa * 0.5)) + self.Ctoned_sa - (self.Cshear_sa * shear_sa)  # - (Cmeta_sa * meta_sa)
        Act_la = 1 / (1 + np.exp(-Stone_la))  # The activation of the smooth muscles....
        Act_sa = 1 / (1 + np.exp(-Stone_sa))
        p1_la = (lam_D_la - 1) * self.Cpassd_la
        p1_sa = (lam_D_sa - 1) * self.Cpassd_sa
        p2_la = np.exp(p1_la)
        p2_sa = np.exp(p1_sa)
        Tpass_la = self.Cpass_la * p2_la
        Tpass_sa = self.Cpass_sa * p2_sa
        s1_la = (lam_D_la - self.Cactd_la) / (self.Cactdd_la)
        s1_sa = (lam_D_sa - self.Cactd_sa) / (self.Cactdd_sa)
        s2_la = math.pow(s1_la, 2)
        s2_sa = math.pow(s1_sa, 2)
        s3_la = np.exp(-s2_la)
        s3_sa = np.exp(-s2_sa)
        Tact_la = self.Cact_la * s3_la  # The activation of the active tension...
        Tact_sa = self.Cact_sa * s3_sa
        Ttot_la = Tpass_la + (Tact_la * Act_la)
        Ttot_sa = Tpass_sa + (Tact_sa * Act_sa)
        return (Ttot_la) - (Pmid_la * Diam_la * 0.5), (Ttot_sa) - (Pmid_sa * Diam_sa * 0.5)

    # This function returns the saturation of Oxygen as a function of length....
    def Saturation(self,x, *params1):
        D1, D2, Q_tot = params1
        Pw2 = (D1, D2, Q_tot)
        if (0 < x <= 0.61):
            return 0.9731
        elif (0.61 < x <= 1.2):
            D = D1
            S_i = 0.9731
            Q = Q_tot
            X_i = 0.61 * 1e-2
        elif (1.2 < x <= 1.56):
            D = D2
            S_i = self.Saturation(1.2, *Pw2)
            Q = Q_tot / 143
            X_i = 1.2 * 1e-2
        elif (1.56 < x <= 1.61):
            D = 6 * 1e-6
            S_i = self.Saturation(1.56, *Pw2)
            Q = Q_tot / 5515
            X_i = 1.56 * 1e-2
        elif (1.61 < x <= 1.97):
            return self.Saturation(1.61, *Pw2)
        else:
            return self.Saturation(1.97, *Pw2)
        Aw1 = np.pi * ((self.n_a * self.l_a) + (self.n_la * self.l_la) + (self.n_sa * self.l_sa) + (self.n_c * self.l_c))
        Aw2 = 2 * np.pi * ((self.n_a * self.l_a * self.Diam_a * 0.5) + (self.n_la * self.l_la * D1 * 0.5) + (self.n_sa * self.l_sa * D2 * 0.5) + (self.n_c * self.l_c * self.Diam_c * 0.5))
        Aw3 = np.pi * ((self.n_a * self.l_a * self.Diam_a * self.Diam_a * 0.25) + (self.n_la * self.l_la * D1 * D1 * 0.25) + (self.n_sa * self.l_sa * D2 * D2 * 0.25) + (self.n_c * self.l_c * self.Diam_c * self.Diam_c * 0.25))
        Aw3 = Aw3 + np.pi * ((self.n_sv * self.l_sv * self.Diam_sv * self.Diam_sv * 0.25) + (self.n_lv * self.l_lv * self.Diam_lv * self.Diam_lv * 0.25) + (self.n_v * self.l_v * self.Diam_v * self.Diam_v * 0.25))
        Aw3 = Aw3 - ((self.n_c * self.l_c) / (500 * 1e06))
        depth = (-Aw2 + np.sqrt((Aw2 * Aw2) - (4 * Aw1 * Aw3))) / (2 * Aw1)
        q1 = 0.25 * np.pi * self.M_0 * (((D + (2 * depth)) * (D + (2 * depth))) - (D * D))
        x = x * 1e-2
        A1 = S_i - (q1 / (Q * self.c0 * self.H_D)) * (x - X_i)
        if (A1 <= 0):
            return 0  # The Saturation is being returned based on the equation given in the report...
        else:
            return A1

    # This calculates the consumption of ATP as a function of x(distance along the length of the vessel)
    def Consumption(self, x, *params2):
        D1, D2, Q_tot = params2
        Pw = (D1, D2, Q_tot)
        X = x * 1e-2
        if (0 < x <= 0.61):
            C_i = 0.5051 * 1e-3
            return C_i
        elif (0.61 < x <= 1.2):
            D = D1
            C_i = 0.5 * 1e-3
            S_i = 0.97
            Q = Q_tot
            X_i = 0.61 * 1e-2
        elif (1.2 < x <= 1.56):
            D = D2
            C_i = self.Consumption(1.2, *Pw)
            S_i = self.Saturation(1.2, *Pw)
            Q = Q_tot / 143
            X_i = 1.2 * 1e-2
        elif (1.56 < x <= 1.61):
            D = 6 * 1e-6
            C_i = self.Consumption(1.56, *Pw)
            S_i = self.Saturation(1.56, *Pw)
            Q = Q_tot / 5515
            X_i = 1.56 * 1e-2
        elif (1.61 < x <= 1.97):
            D = 23.5 * 1e-6
            C_i = self.Consumption(1.61, *Pw)
            S_i = self.Saturation(1.61, *Pw)
            Q = Q_tot / 143
            X_i = 1.61 * 1e-2

        else:
            D = 119.1 * 1e-6
            C_i = self.Consumption(1.97, *Pw)
            S_i = self.Saturation(1.61, *Pw)
            Q = Q_tot
            X_i = 1.97 * 1e-2

        if (self.Saturation(x, *Pw) <= 0):
            gamma1 = (self.k_d * (np.pi) * D) / (0.6 * Q)
            alpha1 = ((self.H_T * self.R_0) / (4 * self.k_d)) * ((D * (1 - (self.R_1 * S_i))))
            A1 = (C_i - alpha1) / (np.exp(-gamma1 * X_i))
            C = A1 * np.exp(-gamma1 * X) + alpha1
            return C
        else:
            Aw1 = np.pi * ((self.n_a * self.l_a) + (self.n_la * self.l_la) + (self.n_sa * self.l_sa) + (self.n_c * self.l_c))
            Aw2 = 2 * np.pi * ((self.n_a * self.l_a * self.Diam_a * 0.5) + (self.n_la * self.l_la * D1 * 0.5) + (self.n_sa * self.l_sa * D2 * 0.5) + (self.n_c * self.l_c * self.Diam_c * 0.5))
            Aw3 = np.pi * ((self.n_a * self.l_a * self.Diam_a * self.Diam_a * 0.25) + (self.n_la * self.l_la * D1 * D1 * 0.25) + (self.n_sa * self.l_sa * D2 * D2 * 0.25) + (self.n_c * self.l_c * self.Diam_c * self.Diam_c * 0.25))
            Aw3 = Aw3 + np.pi * ((self.n_sv * self.l_sv * self.Diam_sv * self.Diam_sv * 0.25) + (self.n_lv * self.l_lv * self.Diam_lv * self.Diam_lv * 0.25) + (self.n_v * self.l_v * self.Diam_v * self.Diam_v * 0.25))
            Aw3 = Aw3 - ((self.n_c * self.l_c) / (500 * 1e06))
            depth = (-Aw2 + np.sqrt((Aw2 * Aw2) - (4 * Aw1 * Aw3))) / (2 * Aw1)
            q1 = 0.25 * np.pi * self.M_0 * (((D + (2 * depth)) * (D + (2 * depth))) - (D * D))
            alpha = ((self.H_T * self.R_0) / (4 * self.k_d)) * ((D * (1 - (self.R_1 * S_i))) - (((1 - self.H_D) * self.R_1 * q1) / ((np.pi) * self.c0 * self.H_D * self.k_d)))
            beta = (D * self.H_T * self.R_0 * self.R_1 * q1) / (4 * Q * self.c0 * self.H_D * self.k_d)
            gamma = (self.k_d * (np.pi) * D) / (0.6 * Q)
            C = alpha + (beta * (X - X_i)) + np.exp(gamma * (X_i - X)) * (C_i - alpha)
            return C

    # Function inside integral
    def Consumption2(self,x, *params3):
        x11, D1, D2, Q_tot = params3
        params2 = (D1, D2, Q_tot)
        C = np.exp(-(x - x11) / (self.L_0)) * self.Consumption(x * 100, *params2)
        return C

    # This function calculates the SCR value for each segment...
    def SCR(self,x, *params4):
        x1 = x * 1e-2
        params41 = (x1, *params4)
        SCR1, error = integrate.quad(self.Consumption2, x1, self.xend_lv, args=params41)
        # print('error is = ',error)
        return SCR1



    def Calculate(self):

        Diam_sa = 30 * 1e-6  # Diameter of the small arteriole which is given as initial guess to the fsolve routine...
        Diam_la = 100 * 1e-6  # Diameter of the large arteriole which is given as initial guess to the fsolve routine...
        # Assumed to be 14 mmHg.
        D = [Diam_la, Diam_sa]
        Pres = 45  ## For the passive myogenic shear and metabolic case...
        Tlac = 100 * 133.33 * self.Diam_lac * 0.5
        meta_la = ((self.Cmyo_la * Tlac) - (self.Cshear_la * 5.5) + self.Ctonedd_la) / self.Cmeta_la
        print(meta_la)
        Tsac = 100 * 133.33 * self.Diam_sac * 0.5
        meta_sa = ((self.Cmyo_sa * Tsac) - (self.Cshear_sa * 5.5) + self.Ctonedd_sa) / self.Cmeta_sa
        print(meta_sa)
        while (Pres <= 200):
            # List of parameters of the large arteriole...
            # Solve for the diameter of the large arteriole...
            D = fsolve(self.Tension2, D,args=Pres)
            k1, k2 = self.Tension2(D,Pres)
            # print('The value of error for large arteriole = ',abs(k1),' and the value is ',D[0],'Pressure is = ',Pres)
            # print('The value of error for small arteriole = ',abs(k2),'and the value is ',D[1],'Pressure is = ',Pres)
            # Calculate the optimum diameter to reduce the number of fluctuations...
            if (abs(k1) > 1e-08 or abs(k2) > 1e-07):
                n = 1
                while (True):
                    print('Hitting a hole at a diameter (large) = ', D[0], ' Pressure  = ', Pres)
                    D1 = (self.Diameter_la[-1] * 1e-6) + (n * 1e-6)
                    D2 = (self.Diameter_sa[-1] * 1e-6) + (n * 1e-6)
                    D = [D1, D2]
                    D = fsolve(self.Tension2, D,args=Pres)
                    k4, k5 = self.Tension2(D)
                    print('error is = ', k4, ' and ', k5)
                    if (abs(k4) < 1e-4 and abs(k5) < 1e-4):
                        break
                    n += 1

            # Store the pressures and diameters....
            self.Pressure_in.append(Pres)
            self.Diameter_la.append(D[0] * 1e06)
            self.Diameter_sa.append(D[1] * 1e06)
            Pres += 1

        k = 0
        # Loop to calculate the perfusion...
        Pressure_la = []
        Pressure_sa = []
        while (k < len(self.Diameter_la)):
            gradP_tot = (self.Pressure_in[k] - 12.91) * 133.33  # Pressure converted from mmHg to N/m^2 by multiplying with 133.33
            Resistance_la = self.compartment_resistance(self.vis_la, self.l_la, self.Diameter_la[k] * 1e-6, self.n_la)
            Resistance_sa = self.compartment_resistance(self.vis_sa, self.l_sa, self.Diameter_sa[k] * 1e-6, self.n_sa)
            Resistance_c = self.compartment_resistance(self.vis_c, self.l_c, self.Diam_c, self.n_c)
            Resistance_lv = self.compartment_resistance(self.vis_lv, self.l_lv, self.Diam_lv, self.n_lv)
            Resistance_sv = self.compartment_resistance(self.vis_sv, self.l_sv, self.Diam_sv, self.n_sv)
            Resistance_a = self.compartment_resistance(self.vis_a, self.l_a, self.Diam_a, self.n_a)
            Resistance_v = self.compartment_resistance(self.vis_v, self.l_v, self.Diam_v, self.n_v)
            Resistance_total = Resistance_sa + Resistance_la + Resistance_c + Resistance_lv + Resistance_sv
            Resistance_total = Resistance_total + Resistance_a  # Calculate the total resistance....
            # print(Resistance_total)
            Q_tot = gradP_tot / Resistance_total
            Diam_la = self.Diameter_la[k] * 1e-6
            Diam_sa = self.Diameter_sa[k] * 1e-6
            P1 = (self.Pressure_in[k] * 133.333) - (Q_tot * Resistance_a)
            P2 = (self.Pressure_in[k] * 133.333) - (Q_tot * Resistance_a) - (Q_tot * Resistance_la)
            Pmid_la = (P1 + P2) * 0.5
            P22 = (self.Pressure_in[k] * 133.333) - (Q_tot * Resistance_a) - (Q_tot * Resistance_la) - (Q_tot * Resistance_sa)
            P11 = (self.Pressure_in[k] * 133.333) - (Q_tot * Resistance_la) - (Q_tot * Resistance_a)
            Pmid_sa = (P11 + P22) * 0.5
            self.perfusion.append(Q_tot)  # The perfusion values are to be collected...
            Pressure_la.append(Pmid_la / 133.33)
            Pressure_sa.append(Pmid_sa / 133.33)
            if (self.Pressure_in[k] == 100):
                perfuse_100 = 6.939226301150093e-11  # The value of this perfusion is around 0.0117 it is quite same with 0.0118.
                # This value is taken from the paper at the control state...
                print('Perfusion=', perfuse_100)

            k += 1

        # Calculate the normalised perfusion...
        self.perfusion_norm = [(k / perfuse_100) for k in self.perfusion]
        return {'Pressure': self.Pressure_in, 'Diameter(LA)': self.Diameter_la, 'Diameter(SA)': self.Diameter_sa, 'Normalised Perfusion': self.perfusion_norm}
    """ 
    Test_Pressure = []
    Test_Pressure1 = []
    Test_Pressure2 = []
    Test_Diameter = []
    Test_Diameter1 = []
    Test_perfusion = []
    # gradP_100 = ((70.9 / 6000) * (Resistance_total100) * vol_100) / 133
    # print(Diameter_sa)
    # print(S3)

    for d in csv.DictReader(open('./perfusion(myo+shear).csv')):
        Test_Pressure1.append(float(d['Pressure']))
        Test_perfusion.append(float(d['Perfusion']))

    # print('Pressure = ', Test_Pressure)
    # print('Diameter = ', Test_Diameter)
    plt.figure(figsize=(10, 9))
    plt.subplot(312)
    plt.xlabel('Pressure(mmHg)')
    plt.ylabel('Perfusion')
    plt.xlim(0, 200)
    plt.ylim(0, 3)
    plt.plot(Test_Pressure1, Test_perfusion, 'b')
    plt.plot(Pressure_in, perfusion_norm, 'r')
    plt.show()
    interpolated1 = np.interp(Pressure_in, Test_Pressure1, Test_perfusion).tolist()
    k = 0
    sum221 = 0.0
    while (k < len(interpolated1)):
        sum221 = sum221 + (abs(interpolated1[k] - perfusion_norm[k]) / (interpolated1[k]))
        if (perfusion_norm[k] > max(Test_perfusion)):
            break
        k += 1
    avg = sum221 / (len(interpolated1))
    print('the error in Perfusion  = ', avg * 100)
    """

