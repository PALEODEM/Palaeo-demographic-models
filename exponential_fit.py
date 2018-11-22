# -*- coding: utf-8 -*-
"""
    File name: calibrate_single_date_mixed_curve.py
    Author: Mario Gutiérrez-Roig
    Date created: 14/07/2016
    Python Version: 3.4.3
    Description: Makes an exponential fit of the SPD with the confidence
    interval. An exponential function is fitted for the whole range, but also
    for other sub-periods (Magdalenian, Younger Dryas and Mesolithic). The 
    output is a plot with the curves and screen text with parameter values.
"""

import numpy as np
from scipy import stats
import pickle

# Input File (Confidence Interval)
input_file = "Simulations/SPD_ConfidenceIntervals_Iberia.pkl"

# Output file plots
output_file = "Plots/SPD_Iberia_exponential_fit.pdf"

# Import Confidence Interval File
with open(input_file, 'rb') as f:
    u = pickle._Unpickler(f)
    u.encoding = 'latin1'
    ci = u.load()

# Define the error of the SPD as the 68.2% (1-sigma) interval
err = 0.5 * (np.array(ci[4])-np.array(ci[1]))

# Import python package for fit 
from scipy.optimize import curve_fit

# Define the function we want to fit
def func(x, lam, C):
    ''' No-normalized exponential function '''
    return C*np.exp(-lam*x)


### EXPONENTIAL FIT WHOLE RANGE ###

# Define the range of the fit
rangemin0 = 8000
rangemax0 = 16600
# Determine the indexes in the array of range bounds
idxmax0 = min(enumerate(ci[0]), key=lambda x: abs(x[1]-rangemax0))[0]
idxmin0 = min(enumerate(ci[0]), key=lambda x: abs(x[1]-rangemin0))[0]

# Perform the fit
popt, pcov = curve_fit(func, ci[0][idxmin0:idxmax0], ci[3][idxmin0:idxmax0], p0=(1e-4, 1))
# Generate the curve resulting from the fit
expfit0 = func(ci[0][idxmin0:idxmax0], *popt)

# Calculate the correlation coefficient and chi-square values
r = np.corrcoef(ci[3][idxmin0:idxmax0],expfit0)[0,1]
chi2 = stats.chisquare(ci[3][idxmin0:idxmax0],expfit0)

# Get the lamda parameter in the exponential and its error
lam = popt[0]
errlam = pcov[0,0]**0.5

# Calculate the % Annual Growth Rate and its error
gr = (np.exp(lam)-1)*100
errgr = lam*np.exp(lam)*errlam*100

# Print the results
print("\n----- EXPONENTIAL FITTING WHOLE RANGE -----")
print(" f(x) = C * exp(-lam * x)")
print(" from ", rangemax0, "BP to", rangemin0, "BP")
print(" lam=",popt[0],"+/-",pcov[0,0]**0.5)
print(" C=",popt[1],"+/-",pcov[1,1]**0.5)
print(" Correlation Coeficient r=",r)
print(" % Annual Growth Rate =",gr,"+/-",errgr)

### EXPONENTIAL FIT MAGDALENIAN ###

# Define the range of the fit
rangemin1 = 12800
rangemax1 = 14860
# Determine the indexes in the array of range bounds
idxmax1 = min(enumerate(ci[0]), key=lambda x: abs(x[1]-rangemax1))[0]
idxmin1 = min(enumerate(ci[0]), key=lambda x: abs(x[1]-rangemin1))[0]

# Perform the fit
popt, pcov = curve_fit(func, ci[0][idxmin1:idxmax1], ci[3][idxmin1:idxmax1], p0=(1e-4, 1), sigma=err[idxmin1:idxmax1])
# Generate the curve resulting from the fit
expfit1 = func(ci[0][idxmin1:idxmax1], *popt)

# Calculate the correlation coefficient and chi-square values
r = np.corrcoef(ci[3][idxmin1:idxmax1], expfit1)[0,1]
chi2 = stats.chisquare(ci[3][idxmin1:idxmax1],expfit1)

# Get the lamda parameter in the exponential and its error
lam = popt[0]
errlam = pcov[0,0]**0.5

# Calculate the % Annual Growth Rate and its error
gr = (np.exp(lam)-1)*100
errgr = lam*np.exp(lam)*errlam*100

# Print the results
print("\n----- EXPONENTIAL FIT MAGDALENIAN -----")
print(" f(x) = C * exp(-lam * x)")
print(" from ", rangemax1, "BP to", rangemin1, "BP")
print(" lam=",popt[0],"+/-",pcov[0,0]**0.5)
print(" C=",popt[1],"+/-",pcov[1,1]**0.5)
print(" Correlation Coeficient r=",r)
print(" % Annual Growth Rate =",gr,"+/-",errgr)

### EXPONENTIAL FIT YOUNGER DRYAS ###

# Define the range of the fit
rangeminYD = 10200
rangemaxYD = 12750
# Determine the indexes in the array of range bounds
idxmaxYD = min(enumerate(ci[0]), key=lambda x: abs(x[1]-rangemaxYD))[0]
idxminYD = min(enumerate(ci[0]), key=lambda x: abs(x[1]-rangeminYD))[0]

# Perform the fit
popt, pcov = curve_fit(func, ci[0][idxminYD:idxmaxYD], ci[3][idxminYD:idxmaxYD], p0=(-1e-3, 1e-4), sigma=err[idxminYD:idxmaxYD])
# Generate the curve resulting from the fit
expfitYD = func(ci[0][idxminYD:idxmaxYD], *popt)

# Calculate the correlation coefficient and chi-square values
r = np.corrcoef(ci[3][idxminYD:idxmaxYD], expfitYD)[0,1]
chi2 = stats.chisquare(ci[3][idxminYD:idxmaxYD],expfitYD)

# Get the lamda parameter in the exponential and its error
lam = popt[0]
errlam = pcov[0,0]**0.5

# Calculate the % Annual Growth Rate and its error
gr = (np.exp(lam)-1)*100
errgr = lam*np.exp(lam)*errlam*100

# Print the results
print("\n----- EXPONENTIAL FIT YOUNGER DRYAS -----")
print(" f(x) = C * exp(-lam * x)")
print(" from ", rangemaxYD, "BP to", rangeminYD, "BP")
print(" lam=",popt[0],"+/-",pcov[0,0]**0.5)
print(" C=",popt[1],"+/-",pcov[1,1]**0.5)
print(" Correlation Coeficient r=",r)
print(" % Annual Growth Rate =",gr,"+/-",errgr)

### EXPONENTIAL FIT MESOLITHIC ###

# Define the range of the fit
rangemin2 = 8000
rangemax2 = 10200
# Determine the indexes in the array of range bounds
idxmax2 = min(enumerate(ci[0]), key=lambda x: abs(x[1]-rangemax2))[0]
idxmin2 = min(enumerate(ci[0]), key=lambda x: abs(x[1]-rangemin2))[0]

# Perform the fit
popt, pcov = curve_fit(func, ci[0][idxmin2:idxmax2], ci[3][idxmin2:idxmax2], p0=(1e-4, 1), sigma=err[idxmin2:idxmax2])
# Generate the curve resulting from the fit
expfit2 = func(ci[0][idxmin2:idxmax2], *popt)

# Calculate the correlation coefficient and chi-square values
r = np.corrcoef(ci[3][idxmin2:idxmax2],expfit2)[0,1]
chi2 = stats.chi2_contingency(ci[3][idxmin2:idxmax2],expfit2)

# Get the lamda parameter in the exponential and its error
lam = popt[0]
errlam = pcov[0,0]**0.5

# Calculate the % Annual Growth Rate and its error
gr = (np.exp(lam)-1)*100
errgr = lam*np.exp(lam)*errlam*100

# Print the results
print("\n----- EXPONENTIAL FIT MESOLITHIC -----")
print(" f(x) = C * exp(-lam * x)")
print(" from ", rangemax2, "BP to", rangemin2, "BP")
print(" lam=",popt[0],"+/-",pcov[0,0]**0.5)
print(" C=",popt[1],"+/-",pcov[1,1]**0.5)
print(" Correlation Coeficient r=",r)
print(" % Annual Growth Rate =",gr,"+/-",errgr)


maxi = 18000
mini = 7000

# Plotting the sum of probilities
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines

fig, ax = plt.subplots(figsize=(12,6))

plt.tick_params(labelsize=14)
plt.xlim(maxi,mini)
#plt.ylim(0,0.0005)
plt.xlabel("Calendar Age (BP)", fontsize=18)
plt.ylabel("Probability", fontsize=18)
plt.fill_between(ci[0], ci[5], ci[1], facecolor='grey', alpha=0.25, edgecolor="none")
plt.fill_between(ci[0], ci[4], ci[2], facecolor='grey', alpha=0.60, edgecolor="none")
plt.plot(ci[0], ci[3], '-',  linewidth=2, color='black')

plt.plot(ci[0][idxmin0:idxmax0], expfit0, '-',  linewidth=3, color='blue')
plt.plot(ci[0][idxmin1:idxmax1], expfit1, '-',  linewidth=3, color='green')
plt.plot(ci[0][idxminYD:idxmaxYD], expfitYD, '-',  linewidth=3, color='cyan')
plt.plot(ci[0][idxmin2:idxmax2], expfit2, '-',  linewidth=3, color='red')

# Uncomment this for logscale
#plt.ylim(1e-5,5e-4)
#plt.yscale('log')

boot = mlines.Line2D([], [], linewidth=1, color='black')
lci95 = mpatches.Patch(color='grey', alpha=0.25)
lci68 = mpatches.Patch(color='grey', alpha=0.60)
exp0 = mlines.Line2D([], [], ls='-', linewidth=3, color='blue')
exp1 = mlines.Line2D([], [], ls='-', linewidth=3, color='green')
expYD = mlines.Line2D([], [], ls='-', linewidth=3, color='cyan')
exp2 = mlines.Line2D([], [], ls='-', linewidth=3, color='red')

plt.legend([boot, lci95, lci68, exp0, exp1, expYD, exp2], ["Population Proxy", "Pop. Prox. 95.4% C.I.", "Pop. Prox. 68.2% C.I.", "Exponential Fit Whole Range", "Exponential Fit Magdalenian", "Exponential Fit Younger Dryas", "Exponential Fit Mesolithic"], loc=2)
plt.savefig(output_file, format='pdf')

