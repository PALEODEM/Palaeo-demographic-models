# -*- coding: utf-8 -*-
"""
    File name: frequencies_histogram.py
    Author: Mario Gutiérrez-Roig
    Date created: 14/07/2016
    Python Version: 3.4.3
    Description: This code plots the frequency histogram of calibrated dates
    medians. Radiocarbon dates are imported from databases. It only takes into
    account filtered dates according an archaeological criteria. Dates with a
    certain percentage of marine source carbon-14 are calibrated with a mixture
    of atmospherical-marine curve.
"""

import numpy as np
import iosacal
import pkg_resources

# Input data file name 

input_file="./Data/All_dataset_F2.csv"

# Output file name

output_file="./Plots/Frequencies_F2_number of dates.pdf"

# Raw data taking some colums of the input file (fields should be delimited by ";")
raw = np.genfromtxt(input_file, delimiter=";", names=True, usecols=("ID_date", "C14Age", "Std", "Type_site", "Marine", "Marine_err", "Reservoir", "Reservoir_err", "Filtered"), dtype=('|S16', float, float, '|S16', float, float, float, float, int), skip_footer=0)

# Consider only the filtered dates 
filtered = [[x[0], x[1], x[2], x[3]] for x in raw if x[4] == 0 and x[8] == 1]
marinedates = [[x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7]] for x in raw if x[4] > 0 and x[8] == 1]

# Radiocarbon dates in iosacal format
r_dates = [iosacal.core.RadiocarbonDetermination(x[1], x[2], x[0]) for x in filtered ]

# Calibrating non marine dates 
caldates = [x.calibrate('intcal13') for x in r_dates]

# Marine dates
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
        calmarine.append(cal_r)

# Define bins
binwidth = 200
bins = np.arange(7000, 20001, binwidth)

#Generating Frequency Histrogram for each site type
freqs = iosacal.core.FreqHist(caldates, bins)


# x-range limits
maxx = 18000
minx = 7500

# Plot Frequency Histogram
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

minorLocator = AutoMinorLocator()

fig, ax = plt.subplots(figsize=(12,6))
plt.xlim(maxx, minx)
plt.tick_params(labelsize=14)
plt.xlabel("Calendar Age (BP)", fontsize=14)
plt.ylabel("Frequency", fontsize=14)
ax.xaxis.set_minor_locator(minorLocator)

# Uncomment in case of log-scale plot
#plt.yscale('log')

# Uncomment this for lines
#plt.plot(freqs.T[0], freqs.T[1], '-',  linewidth=3, color='red')
# Uncomment this for bars
plt.bar(freqs.T[0] - 0.5*binwidth, freqs.T[1], width=binwidth, color='blue', alpha=0.6)

plt.legend(loc=2)
plt.savefig(output_file, format='pdf')

