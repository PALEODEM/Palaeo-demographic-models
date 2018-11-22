# -*- coding: utf-8 -*-
"""
    File name: exponential_null_model_pvalue.py
    Authors: Mario Gutiérrez-Roig and Fabio Silva
    Date created: 25/10/2017
    Python Version: 3.4.3
    Description: Calculates the significance test overall p-value
"""

import numpy as np
import iosacal
import pickle

# Input data file name 
#input_file_0="./Simulations/SPD_Cantabrian_ConfidenceIntervals.pkl"
#input_file_0="./Simulations/SPD_EbroValley_ConfidenceIntervals.pkl"
#input_file_0="./Simulations/SPD_Mediterranean_ConfidenceIntervals.pkl"
#input_file_0="./Simulations/SPD_Interior_ConfidenceIntervals.pkl"
#input_file_0="./Simulations/SPD_Atlantic_ConfidenceIntervals.pkl"
input_file_0="./Simulations/SPD_Iberia_ConfidenceIntervals.pkl"

# Input data file name
input_file_1="./Simulations/Iberia_exponential_bootstrap.pkl"


# Output data file name
#output_file="./Plots/Cantabrian_exponential_nullmodel.pdf"
#output_file="./Plots/EbroValley_exponential_nullmodel.pdf"
#output_file="./Plots/Mediterranean_exponential_nullmodel.pdf"
#output_file="./Plots/All_exponential_nullmodel_prova.pdf"
#output_file="./Plots/All_exponential_nullmodel_pvalue.pdf"
output_file="./Plots/Iberia_exponential_nullmodel_pvalue.pdf"


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

maxi = 22000
mini = 7500


# Calculating p-value
minp = int(max(min(ci[0]), min(ci_exp[0])))
maxp = int(min(max(ci[0]), max(ci_exp[0])))
xx = range(minp, maxp)
ind1 = np.arange(ci[0].shape[0])[np.in1d(ci[0], xx)]
ind2 = np.arange(ci_exp[0].shape[0])[np.in1d(ci_exp[0], xx)]

# Observed statistic for p-value estimation
area_above = ci[3, ind1] - ci_exp[5, ind2]
area_above[ area_above < 0 ] = 0
area_below = ci_exp[1, ind2] - ci[3, ind1]
area_below[ area_below < 0 ] = 0
obs_stat = sum(area_above) + sum(area_below)

# Expected statistic from all the nullmodel bootstraps
mini = np.int(np.min([np.min(x.T[0]) for x in SPDexpboot]))
maxi = np.int(np.max([np.max(x.T[0]) for x in SPDexpboot]))
Nboot = len(SPDexpboot)
yrange=[]
for spd in SPDexpboot:
    y = np.array([ [val[0], val[1]] for val in spd if (val[0] >= mini) and (val[0] <= maxi)])
    yy = np.lib.pad(y.T[1], (np.int(y.T[0][0] - mini), np.int(maxi - y.T[0][-1])), 'constant', constant_values=0)
    yrange.append(yy)
yrange = np.array(yrange).T

nsims = SPDexpboot.size
np_sims = 0
for i in range(1, nsims):
  area_above = yrange[:,i] - ci_exp[5]
  area_above[ area_above < 0 ] = 0
  area_below = ci_exp[1] - yrange[:,i]
  area_below[ area_below < 0 ] = 0
  exp_stat = sum(area_above) + sum(area_below)
  if exp_stat >= obs_stat:
      np_sims = np_sims + 1

# Estimated p-value
p.val = (np_sims + 1 ) / ( nsims + 1 )

print("p-value = ", p.val)

