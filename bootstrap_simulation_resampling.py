# -*- coding: utf-8 -*-
"""
    File name: bootstrap_simulation_resampling.py
    Author: Mario Gutiérrez-Roig
    Date created: 14/07/2016
    Python Version: 3.4.3
    Description: This code callibrates and generates a SPD from dates stored
    in datasets. Takes into account filtered dates, marine source carbon, and
    taphonomic effects (Surovell's curve) depending on site type (Open,
    Shelter, Cave). After calculating the SPD, a bootstrap resampling
    simulations are performed taking into account the original SPD. The
    collection of simulated SPDs is stored in pickle format as an output.
"""

import numpy as np
import iosacal
import pkg_resources
import pickle

# Input data file name 
#input_file="./Data/Cantabrian_dataset_F2.csv"
#input_file="./Data/EbroValley_dataset_F2.csv"
#input_file="./Data/Mediterranean_dataset_F2.csv"
#input_file="./Data/Atlantic_dataset_F2.csv"
#input_file="./Data/Filter2/Interior_dataset_F2.csv"
input_file="./Data/Filter2/All_dataset_F2.csv"

# Output data file name
#output_file="./Simulations/Cantabrian_bootstrap.pkl"
#output_file="./Simulations/EbroValley_bootstrap.pkl"
#output_file="./Simulations/Mediterranean_bootstrap.pkl"
#output_file="./Simulations/Atlantic_bootstrap.pkl"
#output_file="./Simulations/Interior_bootstrap.pkl"
output_file="./Simulations/All_bootstrap.pkl"

# Raw data taking some colums of the input file (fields should be delimited by ";")
raw = np.genfromtxt(input_file, delimiter=";", names=True, usecols=("ID_date", "C14Age", "Std", "Type_site", "Marine", "Marine_err", "Reservoir", "Reservoir_err", "Filtered"), dtype=('|S32', float, float, '|S16', float, float, float, float, int), skip_footer=0)

# Take only the filtered dates 
filtered = [[x[0], x[1], x[2], x[3]] for x in raw if x[4] == 0 and x[8] == 1]
# Take dates with marine curve separatedly
marinedates = [[x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7]] for x in raw if x[4] > 0 and x[8] == 1]
# Generate a list with standard deviations
stds = [x[2] for x in raw if x[8] == 1]

# Radiocarbon dates seggregated by site type
open_dates = [iosacal.core.RadiocarbonDetermination(x[1], x[2], x[0]) for x in filtered if x[3].decode('utf-8') == "Open"]
shelter_dates = [iosacal.core.RadiocarbonDetermination(x[1], x[2], x[0]) for x in filtered if x[3].decode('utf-8') == "Shelter"]
cave_dates = [iosacal.core.RadiocarbonDetermination(x[1], x[2], x[0]) for x in filtered if x[3].decode('utf-8') == "Cave"]

# Calibrate dates
calopen = [x.calibrate('intcal13', norm=True) for x in open_dates]
calshelter = [x.calibrate('intcal13', norm=True) for x in shelter_dates]
calcave = [x.calibrate('intcal13', norm=True) for x in cave_dates]

# Calculate SPDs
SPDopen = iosacal.core.SPD(calopen, norm=False)
SPDshelter = iosacal.core.SPD(calshelter, norm=False)
SPDcave = iosacal.core.SPD(calcave, norm=False)

# Perform the taphonomic correction
SPDopentaph = SPDopen.taphcorr(function='Surovell', factor=1)
SPDsheltertaph = SPDshelter.taphcorr(function='Surovell', factor=0)
SPDcavetaph = SPDcave.taphcorr(function='Surovell', factor=0)

# Perform the same analysis for marine dates
if len(marinedates) == 0:
    # Gets here in case there are no marine dates
    SPDfinal = iosacal.core.spdsum([SPDcavetaph, SPDsheltertaph, SPDopentaph])
    SPDfinal.ndates = SPDcavetaph.ndates + SPDsheltertaph.ndates + SPDopentaph.ndates
else:
    # Calibrate Marine dates
    calmarine = []
    for x in marinedates:
        r = iosacal.core.RadiocarbonDetermination(x[1], x[2], x[0])
        
        # Indicate the calibration curve
        curve_name = "intcal13"
        # Determine the path to the calibration curve file
        curve_path = pkg_resources.resource_filename("iosacal", "data/%s.14c" % curve_name)
        # CalibrationCurve format of iosacal
        calibration_curve = iosacal.core.CalibrationCurve(curve_path)
        # Mix it with other curve
        calibration_curve.mixing("marine13", P=x[4], D=x[5], deltaR=x[6], err_deltaR=x[7])
        
        cal_r = r.calibrate(calibration_curve)
        calmarine.append(cal_r)
    
    SPDmarine = iosacal.core.SPD(calmarine, norm=False)
    
    SPDmarinetaph = SPDmarine.taphcorr(function='Surovell', factor=1)
    
    SPDfinal = iosacal.core.spdsum([SPDcavetaph, SPDsheltertaph, SPDopentaph, SPDmarinetaph])
    SPDfinal.ndates = SPDcavetaph.ndates + SPDsheltertaph.ndates + SPDopentaph.ndates + SPDmarinetaph.ndates

# Number of bootstrap resamplings
Nboot = 1250

SPDboot = []
# Bootsrap loop
for i in np.arange(Nboot):
    print("Calculating... Iteration", i+1, "of", Nboot)    
    # Simulating the SPD drawing stds from stds list    
    simSPD = SPDfinal.simSPD("intcal13", c14std=stds, seed=np.random.randint(1000))
    SPDboot.append(simSPD)

# Save file
with open(output_file, 'wb') as output:
    pickle.dump(SPDboot, output, -1)



     
