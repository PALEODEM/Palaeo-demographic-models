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

input_file="./Data/All_dataset_F2.csv"

# Output file name

output_file="./Plots/SPD_All_taphonomic.pdf"

# Output file name 2 (pickle file)

output_file_2="./Simulations/SPD_All.pkl"

# Raw data taking some columns of the input file (fields should be delimited by ";")
raw = np.genfromtxt(input_file, delimiter=";", names=True, usecols=("ID_date", "C14Age", "Std", "Type_site", "Marine", "Marine_err", "Reservoir", "Reservoir_err", "Filtered"), dtype=('|S16', float, float, '|S16', float, float, float, float, int), skip_footer=0)

# Consider only the filtered dates 
filtered = [[x[0], x[1], x[2], x[3]] for x in raw if x[4] == 0 and x[8] == 1]
marinedates = [[x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7]] for x in raw if x[4] > 0 and x[8] == 1]

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

# Performs the taphonomic correction
SPDopentaph = SPDopen.taphcorr(function='Surovell', factor=1)
SPDsheltertaph = SPDshelter.taphcorr(function='Surovell', factor=0)
SPDcavetaph = SPDcave.taphcorr(function='Surovell', factor=0)

# Performs the same analysis for marine dates
if len(marinedates) == 0:
    # Gets here in case there are no marine dates
    # Sum all SPDs from different site types into one single SPD
    SPDmixedtaph = iosacal.core.spdsum([SPDcavetaph, SPDsheltertaph, SPDopentaph])
    SPDmixedtaph.ndates = SPDcavetaph.ndates + SPDsheltertaph.ndates + SPDopentaph.ndates
else:
    # Calibrate Marine dates
    calmarine = []
    for x in marinedates:
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
        calmarine.append(cal_r)
    
    # Calculates the SPD from the collection of calibrated marine dates
    SPDmarine = iosacal.core.SPD(calmarine, norm=False)
    
    # Taphonomic correction for marine dates
    SPDmarinetaph = SPDmarine.taphcorr(function='Surovell', factor=1)
    
    # Sum all SPDs from different site types into one single SPD
    SPDmixedtaph = iosacal.core.spdsum([SPDcavetaph, SPDsheltertaph, SPDopentaph, SPDmarinetaph])
    SPDmixedtaph.ndates = SPDcavetaph.ndates + SPDsheltertaph.ndates + SPDopentaph.ndates + SPDmarinetaph.ndates

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


plt.fill_between(SPDmixedtaph.T[0], SPDmixedtaph.T[1], facecolor="blue", alpha="0.25", edgecolor="none")
plt.savefig(output_file, format='pdf')
