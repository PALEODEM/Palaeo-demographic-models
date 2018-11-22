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

#Input data file name 
input_file="Data/Iberia.csv"

# Output data file name
output_file="Simulations/SPD_Iberia_bootstrap_1000.pkl"

# Raw data taking some colums of the input file (fields should be delimited by ";")
raw = np.genfromtxt(input_file, delimiter=";", names=True, usecols=("ID_date", "C14Age", "Std", "Type_site", "Marine", "Marine_err", "Reservoir", "Reservoir_err", "Filtered"), dtype=('|S16', float, float, '|S16', float, float, float, float, int), skip_footer=0)

# Consider only the filtered dates 
filtered = [[x[0], x[1], x[2], x[3]] for x in raw if x[4] == 0 and x[8] == 1]
marinedates = [[x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7]] for x in raw if x[4] > 0 and x[8] == 1]

stds = [x[2] for x in raw if x[8] == 1]

# Marine radiocarbon dates seggreated by site type 
open_marine = [[x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7]] for x in raw if (x[4] > 0 and x[8] == 1) and x[3].decode('utf-8') == "Open"]
shelter_marine = [[x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7]] for x in raw if x[4] > 0 and x[8] == 1 and x[3].decode('utf-8') == "Shelter"]
cave_marine = [[x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7]] for x in raw if x[4] > 0 and x[8] == 1 and x[3].decode('utf-8') == "Cave"]

# Radiocarbon dates seggregated by site type
open_dates = [iosacal.core.RadiocarbonDetermination(x[1], x[2], x[0]) for x in filtered if x[3].decode('utf-8') == "Open"]
shelter_dates = [iosacal.core.RadiocarbonDetermination(x[1], x[2], x[0]) for x in filtered if x[3].decode('utf-8') == "Shelter"]
cave_dates = [iosacal.core.RadiocarbonDetermination(x[1], x[2], x[0]) for x in filtered if x[3].decode('utf-8') == "Cave"]

# Calibrate dates
calopen = [x.calibrate('intcal13', norm=True) for x in open_dates]
calshelter = [x.calibrate('intcal13', norm=True) for x in shelter_dates]
calcave = [x.calibrate('intcal13', norm=True) for x in cave_dates]


# Calculate SPDs and perform the taphonomic correction for atmospheric samples
if len(calopen)>0:
    SPDopen = iosacal.core.SPD(calopen, norm=False)
    SPDopentaph = SPDopen.taphcorr(function='Surovell', factor=1)
if len(calshelter)>0:
    SPDshelter = iosacal.core.SPD(calshelter, norm=False)
    SPDsheltertaph = SPDshelter.taphcorr(function='Surovell', factor=0)
if len(calcave)>0:
    SPDcave = iosacal.core.SPD(calcave, norm=False)
    SPDcavetaph = SPDcave.taphcorr(function='Surovell', factor=0)

# Performs the same analysis for marine dates, talking each taphanomic context in turn
if len(open_marine)>0:
    # Calibrate open marine dates
    calopenmarine = []
    for x in open_marine:
        r = iosacal.core.RadiocarbonDetermination(x[1], x[2], x[0])
        
        # Define the main calibration curve
        curve_name = "intcal13"
        # Determine the path to the calibration curve file
        curve_path = pkg_resources.resource_filename("iosacal", "data/%s.14c" % curve_name)
        # CalibrationCurve format of iosacal
        calibration_curve = iosacal.core.CalibrationCurve(curve_path)
        # Mix it with secondary curve with a certain percentage P +/- D and Rerservoir effect deltaR +/- err_deltaR
        calibration_curve.mixing("marine13", P=x[4], D=x[5], deltaR=x[6], err_deltaR=x[7])
        
        cal_r = r.calibrate(calibration_curve)
        calopenmarine.append(cal_r)
    
    # Calculates the SPD from the collection of calibrated marine dates
    SPDopenmarine = iosacal.core.SPD(calopenmarine, norm=False)
    
    # Taphonomic correction for marine dates from open sites
    SPDopenmarinetaph = SPDopenmarine.taphcorr(function='Surovell', factor=1)

if len(shelter_marine)>0:
    # Calibrate Marine dates
    calsheltermarine = []
    for x in shelter_marine:
        r = iosacal.core.RadiocarbonDetermination(x[1], x[2], x[0])
        
        # Define the main calibration curve
        curve_name = "intcal13"
        # Determine the path to the calibration curve file
        curve_path = pkg_resources.resource_filename("iosacal", "data/%s.14c" % curve_name)
        # CalibrationCurve format of iosacal
        calibration_curve = iosacal.core.CalibrationCurve(curve_path)
        # Mix it with secondary curve with a certain percentage P +/- D and Rerservoir effect deltaR +/- err_deltaR
        calibration_curve.mixing("marine13", P=x[4], D=x[5], deltaR=x[6], err_deltaR=x[7])
        
        cal_r = r.calibrate(calibration_curve)
        calsheltermarine.append(cal_r)
    
    # Calculates the SPD from the collection of calibrated marine dates
    SPDsheltermarine = iosacal.core.SPD(calsheltermarine, norm=False)
    
    # Taphonomic correction for marine dates
    SPDsheltermarinetaph = SPDsheltermarine.taphcorr(function='Surovell', factor=0)

if len(cave_marine)>0:
    # Calibrate Marine dates
    calcavemarine = []
    for x in cave_marine:
        r = iosacal.core.RadiocarbonDetermination(x[1], x[2], x[0])
        
        # Define the main calibration curve
        curve_name = "intcal13"
        # Determine the path to the calibration curve file
        curve_path = pkg_resources.resource_filename("iosacal", "data/%s.14c" % curve_name)
        # CalibrationCurve format of iosacal
        calibration_curve = iosacal.core.CalibrationCurve(curve_path)
        # Mix it with secondary curve with a certain percentage P +/- D and Rerservoir effect deltaR +/- err_deltaR
        calibration_curve.mixing("marine13", P=x[4], D=x[5], deltaR=x[6], err_deltaR=x[7])
        
        cal_r = r.calibrate(calibration_curve)
        calcavemarine.append(cal_r)
    
    # Calculates the SPD from the collection of calibrated marine dates
    SPDcavemarine = iosacal.core.SPD(calcavemarine, norm=False)
    
    # Taphonomic correction for marine dates
    SPDcavemarinetaph = SPDcavemarine.taphcorr(function='Surovell', factor=0)


# Sum all SPDs from different site types into one single SPD
# Manually adjust these lines if certain SPDs were not calculated because there were no dates of that type in the input data
# The full list of different spds are as follows:
#   SPDmixedtaph = iosacal.core.spdsum([SPDcavetaph, SPDsheltertaph, SPDopentaph, SPDopenmarinetaph, SPDsheltermarinetaph, SPDcavemarinetaph])
#   SPDmixedtaph.ndates = SPDcavetaph.ndates + SPDsheltertaph.ndates + SPDopentaph.ndates + SPDcavemarinetaph.ndates + SPDsheltermarinetaph.ndates + SPDopenmarinetaph.ndates

SPDmixedtaph = iosacal.core.spdsum([SPDcavetaph, SPDsheltertaph, SPDopentaph, SPDopenmarinetaph, SPDsheltermarinetaph, SPDcavemarinetaph])
SPDmixedtaph.ndates = SPDcavetaph.ndates + SPDsheltertaph.ndates + SPDopentaph.ndates + SPDcavemarinetaph.ndates + SPDsheltermarinetaph.ndates + SPDopenmarinetaph.ndates

# Number of bootstrap resamplings
Nboot = 100

SPDboot = []
# Bootsrap loop
for i in np.arange(Nboot):
    print("Calculating... Iteration", i+1, "of", Nboot)    
    # Simulating the SPD drawing stds from stds list    
    simSPD = SPDmixedtaph.simSPD("intcal13", c14std=stds, seed=np.random.randint(1000))
    SPDboot.append(simSPD)

# Save file
with open(output_file, 'wb') as output:
    pickle.dump(SPDboot, output, -1)

#with open('C:/Users/mgutierrez/Desktop/Simulations/SPDoriginal.pkl', 'wb') as output:
#    pickle.dump(SPDfinal, output, -1)

