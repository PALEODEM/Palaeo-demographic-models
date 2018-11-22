# -*- coding: utf-8 -*-
"""
    File name: SPD_curve_taphonomic_correction.py
    Author: Mario Gutiérrez-Roig
    Date created: 14/07/2016
    Python Version: 3.4.3
    Description: This code generates a plot of the Sum of Probabilities
    Distribution (SPD) of a collection of dates stored in datasets called as
    an input. It onky takes into account filtered dates according an
    archaeological criteria. Dates with a certain percentage of marine source 
    carbon are calibrated with a mixture of atmospherical-marine curve.
    Taphonomical correction is applied by using Surovell's curve and according
    to different origin of the samples (Open sites, Shelters or Caves).
"""

import numpy as np
import iosacal
import pkg_resources
import pickle

# Input data file path 
input_file="Data/Iberia.csv"

# Output file name
output_file="Plots/SPD_Iberia.pdf"

# Output file name 2 (pickle file)
output_file_2="Simulations/SPD_Iberia.pkl"

# Raw data taking some colums of the input file (fields should be delimited by ";")
raw = np.genfromtxt(input_file, delimiter=";", names=True, usecols=("ID_date", "C14Age", "Std", "Type_site", "Marine", "Marine_err", "Reservoir", "Reservoir_err", "Filtered"), dtype=('|S16', float, float, '|S16', float, float, float, float, int), skip_footer=0)

# Consider only the filtered dates 
filtered = [[x[0], x[1], x[2], x[3]] for x in raw if x[4] == 0 and x[8] == 1]
marinedates = [[x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7]] for x in raw if x[4] > 0 and x[8] == 1]

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


# Save final SPD curve into pickle file
with open(output_file_2, 'wb') as output:
    pickle.dump(SPDmixedtaph, output, -1)

# x-range limits
maxx = 18000
minx = 7500

# Plot script #

import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

minorLocator = AutoMinorLocator()

fig, ax = plt.subplots(figsize=(12,6))

plt.tick_params(labelsize=14)
plt.xlim(maxx,minx)
plt.xlabel("Calendar Age (BP)", fontsize=14)
plt.ylabel("Probability", fontsize=14)
ax.xaxis.set_minor_locator(minorLocator)

# Uncomment in case of log-scale plot
#plt.ylim(5e-6,5e-4)
#plt.yscale('log')

#This plots red and green vertical lines at 8240 and 9300
#plt.plot([8240,8240], [0.0006,0], '-', linewidth=4, color='red')
#plt.plot([9300,9300], [0.0006,0], '-', linewidth=4, color='green')

plt.fill_between(SPDmixedtaph.T[0], SPDmixedtaph.T[1], facecolor="blue", alpha="0.25", edgecolor="none")
plt.plot(SPDmixedtaph.T[0], SPDmixedtaph.T[1], '-', linewidth=1, color='black')
plt.savefig(output_file, format='pdf')