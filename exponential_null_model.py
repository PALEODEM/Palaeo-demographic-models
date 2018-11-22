# -*- coding: utf-8 -*-
"""
    File name: exponential_null_model.py
    Author: Mario Gutiérrez-Roig
    Date created: 14/07/2016
    Python Version: 3.4.3
    Description: This code fits an exponential function on the Population Proxy
    curve. Then, performs a simulation resampling by randomly drawing dates
    from the exponential fitted curve distribution. The number of dates are the
    same than in the original SPD and the standard deviations are also taken
    from the original Datasets. 
    
"""

import numpy as np
import iosacal
import pickle
import pkg_resources

# Input data file name 
input_file_0="Data/Iberia.csv"

# Input File (Original SPD)
input_file_1 = "Simulations/SPD_Iberia.pkl"

# Input File (Confidence Interval)
input_file_2 = "Simulations/SPD_ConfidenceIntervals_Iberia.pkl"

# Output data file name (Simulated data in pickle format)
output_file_1="Simulations/Iberia_Exp.pkl"

# Output file name (Plot)
output_file="Plots/Iberia_Exp.pdf"

# Raw data taking some colums of the input file (fields should be delimited by ";")
raw = np.genfromtxt(input_file_0, delimiter=";", names=True, usecols=("ID_date", "C14Age", "Std", "Type_site", "Marine", "Marine_err", "Reservoir", "Reservoir_err", "Filtered"), dtype=('|S32', float, float, '|S16', float, float, float, float, int), skip_footer=0)

# Generate a list with standard deviations
stds = [x[2] for x in raw if x[8] == 1]

# Importing original SPD
with open(input_file_1, 'rb') as f:
    u = pickle._Unpickler(f)
    u.encoding = 'latin1'
    SPD = u.load() 

# Importing Confidence Intervals
with open(input_file_2, 'rb') as f:
    u = pickle._Unpickler(f)
    u.encoding = 'latin1'
    ci = u.load()

from scipy.optimize import curve_fit

# Defining the range of the curve
rangemin = 7500
rangemax = 22000

# Finding the indexes
idxmax = min(enumerate(ci[0]), key=lambda x: abs(x[1]-rangemax))[0]
idxmin = min(enumerate(ci[0]), key=lambda x: abs(x[1]-rangemin))[0]

# Defining the function to fit
def func(x, lam, C):
    ''' No-normalized exponential function '''
    return C*lam*np.exp(-lam*x)

# Performing the fit
popt, pcov = curve_fit(func, ci[0][idxmin:idxmax], ci[3][idxmin:idxmax], p0=(1e-6, 1))
# Generating the curve
expfit = func(ci[0][idxmin:idxmax], *popt)

# Copying the SPD format from the original SPD 
SPDexp = SPD.copy()
# Resize the SPD
SPDexp.resize(len(ci[0][idxmin:idxmax]), 2, refcheck=False)
# Replace the x-axis and y-axis by the fited curve
SPDexp.T[0] = ci[0][idxmin:idxmax]
SPDexp.T[1] = expfit
# Assigning the corresponding number of dates
SPDexp.ndates = len(stds)

# Number of simulation iterations
Nboot = 1000
SPDexpboot = []
for i in np.arange(Nboot):
    print("Calculating... Iteration", i+1, "of", Nboot)    
    simSPD = SPDexp.simSPD("intcal13", c14std=stds, seed=np.random.randint(1000))
    SPDexpboot.append(simSPD)

# Calculating Confidence Intervals for the exponential function
ci_exp = iosacal.core.confint(SPDexpboot)  

maxi = 18000
mini = 7500

# Save Confidence Intervals in pickle file

with open(output_file_1, 'wb') as output:
    pickle.dump(SPDexpboot, output, -1)


# Plotting the sum of probilities
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines

fig, ax = plt.subplots(figsize=(12,6))

plt.tick_params(labelsize=14)
plt.xlim(maxi,mini)
plt.xlabel("Calendar Age (BP)", fontsize=18)
plt.ylabel("Probability", fontsize=18)
#plt.fill_between(ci[0], ci[1], ci[5], facecolor='grey', alpha=0.25, edgecolor="none")
#plt.fill_between(ci[0], ci[2], ci[4], facecolor='grey', alpha=0.60, edgecolor="none")

#plt.plot(ci[0], ci[3], '-',  linewidth=2, color='black')
plt.plot(SPD.T[0],SPD.T[1], '-',  linewidth=2, color='black')

plt.fill_between(ci_exp[0], ci_exp[1], ci_exp[5], facecolor='blue', alpha=0.15, edgecolor="none")
plt.fill_between(ci_exp[0], ci_exp[2], ci_exp[4], facecolor='blue', alpha=0.35, edgecolor="none")
plt.plot(ci_exp[0], ci_exp[3], '-',  linewidth=2, color='blue')

#plt.ylim(5e-6,1e-3)
#plt.yscale('log')

ori = mlines.Line2D([], [], linewidth=2, color='black')
#lci95 = mpatches.Patch(color='grey', alpha=0.25)
#lci68 = mpatches.Patch(color='grey', alpha=0.60)

exp = mlines.Line2D([], [], linewidth=2, color='blue')
exp95 = mpatches.Patch(color='blue', alpha=0.15)
exp68 = mpatches.Patch(color='blue', alpha=0.35)

plt.legend([ori, exp, exp95, exp68], ["Original SPD", "Exponential null model", "Exponential 95.4% C.I.", "Exponential 68.2% C.I."], loc=2)

plt.savefig(output_file, format='pdf')

