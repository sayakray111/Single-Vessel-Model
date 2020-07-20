# Single-Vessel-Model

This file calculates the Pressure and Diameter relations between the Diameter of the large arteriole and the Pressure at the midpoint of the large arteriole. This file can be modified by a simple way to calculate the perfusion on the vessel. The perfusion_norm term relates the normalised perfusion. Thus to run Normalised perfusion please run the perfusion_norm vs Pressure_in on the last line which says to plot. In that case the figures described in the latex document can be replicated. 

Files Required : DAnometa2vals.mat, Carlson(2008)_dataP.csv, perfusion(passive).csv for Carlson_Passive.py...DAnometa2vals.mat, Carlson(2008)_data.csv, perfusion(myo).csv for the Carlson_Passive_myo.py Carlson_2008(myo+meta+shear).csv and perfusion(myo+shear+meta).csv for the Carlson_metabolic.py
Carlson(Passive).py generates the perfusion and Pressure_Diameter curves for the passive case, Carlson_Passive_myo.py generates the perfusion and Pressure_Diameter curves for the myogenic case, The Carlson_metabolic.py generates the perfusion and Pressure_Diameter curves for the metabolic case...
Libraries Required : Scipy Numpy 

The libraries can easily be installed using pip install either in base environment or virtual environment or using anaconda. 
Installing Numpy (Windows): 

In anaconda either in virtual environment or in base environment...............................
conda install -c anaconda numpy 

Installing Numpy (Linux): 
Using pip  - pip install numpy

Instructions: Git clone the files and then run the output. The code will reproduce the curves given in the Carlson(2008) paper - Theoretical model of blood flow autoregulation: roles of myogenic, shear-dependent, and metabolic responses. The code will give the results depending on the set of parameters given in the paper. 
