# -*- coding: utf-8 -*-
"""
    File name: frequencies_histogram_type_sites.py
    Author: Mario Gutiérrez-Roig
    Date created: 14/07/2016
    Python Version: 3.4.3
    Description: This code plots the frequency histogram of calibrated dates
    medians. The Radiocarbon dates are imported from databases. It only takes into
    account filtered dates according an archaeological criteria. Dates with a
    certain percentage of marine source carbon-14 are calibrated with a mixture
    of atmospherical-marine curve. Dates coming from different sites type
    (Open sites, shelters, caves) are distinguished with different colors. 
"""

import numpy as np
import iosacal

# Input data file name 

input_file="./Data/All_dataset_F2.csv"

# Output file name

output_file="./Plots/Frequencies_F2_type_sites.pdf"

# Raw data taking some colums of the input file (fields should be delimited by ";")
raw = np.genfromtxt(input_file, delimiter=";", names=True, usecols=("ID_date", "C14Age", "Std", "Type_site", "Marine", "Marine_err", "Reservoir", "Reservoir_err", "Filtered"), dtype=('|S16', float, float, '|S16', float, float, float, float, int), skip_footer=0)

# Consider only the filtered dates 
filtered = [[x[0], x[1], x[2], x[3]] for x in raw if x[8] == 1]
unfiltered = [[x[0], x[1], x[2], x[3]] for x in raw if x[8] == 0]

# Radiocarbon dates classified by site type
shelter_dates = [iosacal.core.RadiocarbonDetermination(x[1], x[2], x[0]) for x in filtered if x[3].decode('utf-8') == "Shelter" ]
cave_dates = [iosacal.core.RadiocarbonDetermination(x[1], x[2], x[0]) for x in filtered if x[3].decode('utf-8') == "Cave"]
open_dates = [iosacal.core.RadiocarbonDetermination(x[1], x[2], x[0]) for x in filtered if x[3].decode('utf-8') == "Open"]

# Calibrating dates for each site type
calshelter = [x.calibrate('intcal13') for x in shelter_dates]
calcave = [x.calibrate('intcal13') for x in cave_dates]
calopen = [x.calibrate('intcal13') for x in open_dates]

# Define bins
binwidth = 200
bins=np.arange(7000, 20001, binwidth)

# Generating Frequency Histrogram for each site type
freqs_open = iosacal.core.FreqHist(calopen, bins)
freqs_shelter = iosacal.core.FreqHist(calshelter, bins)
freqs_cave = iosacal.core.FreqHist(calcave, bins)

# x-range limits
maxx = 18000
minx = 7500

import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

minorLocator = AutoMinorLocator()

fig, ax = plt.subplots(figsize=(12,6))

plt.xlim(maxx, minx)
plt.tick_params(labelsize=14)
plt.xlabel("Calendar Age (BP)", fontsize=14)
plt.ylabel("Frequency", fontsize=14)
ax.xaxis.set_minor_locator(minorLocator)

p1 = plt.bar(freqs_cave.T[0] - 0.5*binwidth, freqs_cave.T[1], width=binwidth, color='grey', label="Cave")
p2 = plt.bar(freqs_shelter.T[0] - 0.5*binwidth, freqs_shelter.T[1], width=binwidth, color='#b0b2b1', bottom=freqs_cave.T[1], label="Rock Shelters")
p3 = plt.bar(freqs_open.T[0] - 0.5*binwidth, freqs_open.T[1], width=binwidth, color='white', bottom=freqs_shelter.T[1] + freqs_cave.T[1], label="Open Sites")

plt.legend(loc=2)
plt.savefig(output_file, format='pdf') 

