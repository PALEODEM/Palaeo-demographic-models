# Palaeo-demographic-models

This repository contains all the scripts implemented to be applied in the article “Palaeo-demographic modelling demonstrates a population bottleneck during the Pleistocene-Holocene transition in Iberia".

Scripts are classified as follows:

## 1.- Exploratory analysis of data

*frequencies_histogram.py: 
Plots the frequency histogram of calibrated dates medians. Radiocarbon dates are imported from databases. It only takes into account filtered dates according an archaeological criteria. Dates with a certain percentage of marine source carbon-14 are calibrated with a mixture of atmospherical-marine curve.

*Frequencies_histogram_type_sites: 
Plots the frequency histogram of calibrated dates medians. The Radiocarbon dates are imported from databases. It only takes into account filtered dates according an archaeological criteria. Dates with a certain percentage of marine source carbon-14 are calibrated with a mixture of atmospherical-marine curve. Dates coming from different sites type (Open sites, shelters, caves) are distinguished with different colors.

*SPD_curve.py: 
This code generates a plot of the Sum of Probabilities Distribution (SPD) of a collection of dates stored in datasets called as an input. It only takes into account filtered dates according an archaeological criteria. Dates with a certain percentage of marine source carbon are calibrated with a mixture of atmospherical-marine curve.

*SPD_curve_taphonomic_correction.py: 
This code generates a plot of the Sum of Probabilities Distribution (SPD) of a collection of dates stored in datasets called as an input. It onky takes into account filtered dates according an archaeological criteria. Dates with a certain percentage of marine source carbon are calibrated with a mixture of atmospherical-marine curve. Taphonomical correction is applied by using Surovell's curve and according to different origin of the samples (Open sites, Shelters or Caves).


## 2.- Modeling

*bootstrap_simulation_resampling.py: 
This code calibrates and generates a SPD from dates stored in datasets. Takes into account filtered dates, marine source carbon, and taphonomic effects (Surovell's curve) depending on site type (Open, Shelter, Cave). After calculating the SPD, a bootstrap resampling simulations are performed taking into account the original SPD. The collection of simulated SPDs is stored in pickle format as an output.

*SPD_confidence_intervals.py: 
This code computes and plots the 68.2% and 95.4% confidence intervals of the SPD from the simulated SPDs stored in a pickle file.


## 3.- Hypothesis-testing

*exponential_null_model_pvalue.py: 
this code performs a simulation resampling by randomly drawing dates from a exponential distribution. The number of dates are the same then in the original SPD and the standard deviations are also taken from the original data set. In addition,  it calculates the significance test overall p-value.

## 4.- Model selection

*model_selection.R: 
This code does regression fitting for the six models considered over a range of values of the breakpoints and uses the information-theoretic approach of Burnham and Anderson (2002) to identify the model that best-fits the SPD data.

*plot_best_models.R: 
This code takes the results of model_selection.R to plot the regression fits of the best breakpoint values for each model

*growth_rates.R: 
This code computes the mean growth rate and associated confidence envelope for the model that best-fits the SPD, over the entire range of breakpoint values.

