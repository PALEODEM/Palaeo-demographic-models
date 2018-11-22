# -*- coding: utf-8 -*-
"""
    File name: SPD_confidence_intervals.py
    Author: Mario Gutiérrez-Roig
    Date created: 14/07/2016
    Python Version: 3.4.3
    Description: This code computes and plots the 68.2% and 95.4% confidence
    intervals of the SPD from the simulated SPDs stored in a pickle file.
"""

import numpy as np
import iosacal
import pickle

# Input file 1 (Original SPD)
input_file_1 = "Simulations/SPD_Iberia.pkl"

# Input file 2 (collection of simulated SPDs)
input_file_2 = "Simulations/SPD_Iberia_bootstrap.pkl"

# Output file 1 (pickle file with confidence intervals)
output_file_1 = "Simulations/SPD_ConfidenceIntervals_Iberia.pkl"

# Output file 2 (plot)
output_file_2 = "Plots/SPD_ConfidenceIntervals_Iberia.pdf"


with open(input_file_1, 'rb') as f:
    u = pickle._Unpickler(f)
    u.encoding = 'latin1'
    SPD = u.load()

with open(input_file_2, 'rb') as f:
    u = pickle._Unpickler(f)
    u.encoding = 'latin1'
    SPDboot = u.load()

# Calculating Confidence Intervals
ci = iosacal.core.confint(SPDboot)

# Save Confidence Intervals in pickle file
with open(output_file_1, 'wb') as output:
    pickle.dump(ci, output, -1)

# Defining X-range (Cal BP Age range)
minx = 7500
maxx = 18000

# Plotting the sum of probilities
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
from matplotlib.ticker import AutoMinorLocator

minorLocator = AutoMinorLocator()

fig, ax = plt.subplots(figsize=(12,6))

plt.tick_params(labelsize=14)
ax.xaxis.set_minor_locator(minorLocator)
plt.xlim(maxx,minx)
plt.xlabel("Calendar Age (BP)", fontsize=18)
plt.ylabel("Probability", fontsize=18)
plt.fill_between(ci[0], ci[5], ci[1], facecolor='grey', alpha=0.25, edgecolor="none")
plt.fill_between(ci[0], ci[4], ci[2], facecolor='grey', alpha=0.60, edgecolor="none")
plt.plot(ci[0], ci[3], '-',  linewidth=2, color='black')
plt.plot(SPD.T[0], SPD.T[1], '-',  linewidth=1, color='red')

# Uncomment in case of logscale plot
#plt.ylim(5e-6,1e-3)
#plt.yscale('log')

ori = mlines.Line2D([], [], linewidth=1, color='red')
lci95 = mpatches.Patch(color='grey', alpha=0.25)
lci68 = mpatches.Patch(color='grey', alpha=0.60)
med = mlines.Line2D([], [], linewidth=2, color='black')

plt.legend([ori, lci95, lci68, med], ["Original SPD", "Bootstrap 95.4% C.I.", "Bootstrap 68.2% C.I.", "Bootstrap Median"], loc=2)

plt.savefig(output_file_2, format='pdf')
 