# -*- coding: utf-8 -*-
"""
    File name: SPD_curve.py
    Author: Mario Gutiérrez-Roig
    Date created: 14/07/2016
    Python Version: 3.4.3
    Description: This code generates a plot of the Sum of Probabilities
    Distribution (SPD) of a collection of dates stored in datasets called as
    an input. It onky takes into account filtered dates according an
    archaeological criteria. Dates with a certain percentage of marine source 
    carbon are calibrated with a mixture of atmospherical-marine curve.
"""

import numpy as np
import iosacal
import pkg_resources
import pickle

# Input data file path 
input_file="Data/Iberia.csv"

# Output file name
output_file="Plots/SPD_Iberia_no_taph.pdf"

output_file_2="Simulations/SPD_Iberia_no_taph.pkl"

# Raw data taking some colums of the input file (fields should be delimited by ";")
raw = np.genfromtxt(input_file, delimiter=";", names=True, usecols=("ID_date", "C14Age", "Std", "Type_site", "Marine", "Marine_err", "Reservoir", "Reservoir_err", "Filtered"), dtype=('|S16', float, float, '|S16', float, float, float, float, int), skip_footer=0)

# Consider only the filtered dates 
filtered = [[x[0], x[1], x[2], x[3]] for x in raw if x[4] == 0 and x[8] == 1]
marinedates = [[x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7]] for x in raw if x[4] > 0 and x[8] == 1]

# Radiocarbon dates seggregated by site type
radiocarbon_dates = [iosacal.core.RadiocarbonDetermination(x[1], x[2], x[0]) for x in filtered]

# Calibrate dates
cal_dates = [x.calibrate('intcal13', norm=True) for x in radiocarbon_dates]

# Performs the same analysis for marine dates
if len(marinedates) > 0:
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
        cal_dates.append(cal_r)
    
# Calculate SPDs
SPD = iosacal.core.SPD(cal_dates, norm=False)

# x-range limits
maxx = 18000
minx = 7500

# Save final SPD curve into pickle file
with open(output_file_2, 'wb') as output:
    pickle.dump(SPD, output, -1)

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
#plt.plot([8240,8240], [0.0006,0], '-', linewidth=1, color='red')
#plt.plot([9300,9300], [0.0006,0], '-', linewidth=4, color='green')

plt.fill_between(SPD.T[0], SPD.T[1], facecolor="blue", alpha="0.25", edgecolor="none")
plt.plot(SPD.T[0], SPD.T[1], '-', linewidth=1, color='black')
plt.savefig(output_file, format='pdf')