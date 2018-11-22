# -*- coding: utf-8 -*-
"""
    File name: calibrate_single_date.py
    Author: Mario Gutiérrez-Roig
    Date created: 14/07/2016
    Python Version: 3.4.3
    Description: This code calibrates one single date and plots it in an
    oxcal-like way. 
"""

import iosacal
import pkg_resources

# Output file
output_file = "test.png"

# Date is defined by three fields: 1) Id 2) Uncalibrated Age and 3) Standard error
date=["test", 9000, 50]

# Define a iosacal object called RadiocarbonDetermination for uncallibrated dates
r = iosacal.core.RadiocarbonDetermination(date[1], date[2], date[0])

# Calibrate date with specified mixed curve

# Define the main calibration curve
curve_name = "intcal13"
# Determine the path to the calibration curve file
curve_path = pkg_resources.resource_filename("iosacal", "data/%s.14c" % curve_name)
# CalibrationCurve format of iosacal
calibration_curve = iosacal.core.CalibrationCurve(curve_path)
# Mix it with secondary curve with a certain percentage P +/- D and Rerservoir effect deltaR +/- err_deltaR
calibration_curve.mixing("marine13", P=0, D=0, deltaR=(0), err_deltaR=0)
        
cal_r = r.calibrate(calibration_curve, norm=True)


# Plot output of calibrated date in OxCal-like way generated at provided output path
iosacal.plot.single_plot(cal_r,oxcal=True, output=output_file)
