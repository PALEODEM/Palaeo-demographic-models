# -*- coding: utf-8 -*-
"""
    File name: calibrate_single_date_mixed_curve.py
    Author: Mario Gutiérrez-Roig
    Date created: 14/07/2016
    Python Version: 3.4.3
    Description: This code calibrates one single date with a certain percentage
    of marine source carbon with a reservoir effect. The output plot is in the
    oxcal-like style. 
"""

import iosacal
import pkg_resources

# Output file
output_file = "./Plots/single_plot_date_mixed_curve.png"

# Date is defined by three fields: 1) Id 2) Uncalibrated Age and 3) Standard error
date=["0028", 9015, 45]

# Define a iosacal object called RadiocarbonDetermination for uncallibrated dates
r = iosacal.core.RadiocarbonDetermination(date[1], date[2], date[0])

# Indicate the calibration curve
curve_name = "intcal13"
# Determine the path to the calibration curve file
curve_path = pkg_resources.resource_filename("iosacal", "data/%s.14c" % curve_name)
# CalibrationCurve format of iosacal
calibration_curve = iosacal.core.CalibrationCurve(curve_path)
# Mix it with curve marine13 with a percentage of 30% +/- 5% and Reservoir effect of 10 +/- 6
calibration_curve.mixing("marine13", P=0.3, D=0.05, deltaR=10, err_deltaR=6)

# Calibrate date with mixed curve curve
cal_r = r.calibrate(calibration_curve, norm=True)

# Plot output of calibrated date in OxCal-like way generated at provided output path
iosacal.plot.single_plot(cal_r,oxcal=True, output=output_file)

