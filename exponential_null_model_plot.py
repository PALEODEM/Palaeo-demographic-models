# -*- coding: utf-8 -*-
"""
    File name: exponential_null_model_plot.py
    Author: Mario Gutiérrez-Roig
    Date created: 14/07/2016
    Python Version: 3.4.3
    Description: Plots the Confidence Intervals of exponential null model
    together with the Population Proxy.
    
"""

import numpy as np
import iosacal
import pickle

# Input data file name 
input_file_0="Simulations/SPD_ConfidenceIntervals_Iberia.pkl"

# Intput data file name
input_file_1="/home/fabpsilva/Desktop/Mario_things/Simulations_10_2018/Iberia_Exp.pkl"

# Output file name
output_file="/home/fabpsilva/Desktop/Mario_things/Plots_10_2018/Iberia_exponential_nullmodel.pdf"

# Importing boostrap files
with open(input_file_0, 'rb') as f:
    u = pickle._Unpickler(f)
    u.encoding = 'latin1'
    ci = u.load()

# Importing boostrap files
with open(input_file_1, 'rb') as f:
    u = pickle._Unpickler(f)
    u.encoding = 'latin1'
    SPDexpboot = u.load()
    
ci_exp = iosacal.core.confint(SPDexpboot)

maxi = 18000
mini = 7500

# Plotting the sum of probilities
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
from matplotlib.ticker import AutoMinorLocator

minorLocator = AutoMinorLocator()

fig, ax = plt.subplots(figsize=(12,6))

plt.tick_params(labelsize=14)
plt.xlim(maxi,mini)
plt.ylim(0,0.0004)
plt.xlabel("Calendar Age (BP)", fontsize=18)
plt.ylabel("Probability", fontsize=18)
ax.xaxis.set_minor_locator(minorLocator)

plt.fill_between(ci_exp[0], ci_exp[1], ci_exp[5], facecolor='grey', alpha=0.25, edgecolor="none")
plt.fill_between(ci_exp[0], ci_exp[2], ci_exp[4], facecolor='grey', alpha=0.60, edgecolor="none")
plt.plot(ci_exp[0], ci_exp[3], '-',  linewidth=2, color='black')

#plt.fill_between(ci[0], ci[1], ci[5], facecolor='red', alpha=0.15, edgecolor="none")
#plt.fill_between(ci[0], ci[2], ci[4], facecolor='red', alpha=0.35, edgecolor="none")
plt.plot(ci[0], ci[3], '-',  linewidth=2, color='red')

#plt.ylim(5e-6,1e-3)
#plt.yscale('log')

ori = mlines.Line2D([], [], linewidth=2, color='red')
#lci95 = mpatches.Patch(color='red', alpha=0.15)
#lci68 = mpatches.Patch(color='red', alpha=0.35)

exp = mlines.Line2D([], [], linewidth=2, color='black')
exp95 = mpatches.Patch(color='grey', alpha=0.25)
exp68 = mpatches.Patch(color='grey', alpha=0.60)

plt.legend([ori, exp, exp95, exp68], ["Population Proxy", "Exponential null model", "Exponential 95.4% C.I.", "Exponential 68.2% C.I."], loc=2)

plt.savefig(output_file, format='pdf')

