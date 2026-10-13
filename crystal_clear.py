# -*- coding: utf-8 -*-
"""
Created on Fri Jan 16 11:45:08 2026
@author: miles dougan

Crystal clear: A shortcut to peak picking and indexing in powder x-ray diffraction (PXRD)

"""

# Modules
import math
import numpy as np 
import matplotlib.pyplot as plt
import peak_picking_functions as ppf
import argparse as arg

# ================= USER INPUT =================
# Enter the names of your input files in the terminal. 
# These files should be located in the same directory as this script.

parse = arg.ArgumentParser(
    description = "Crystal Clear Program"
    )

# Required arguments (input files)
parse.add_argument("exp_file", help = "Experimental PXRD data file (CSV)")
parse.add_argument("ref_file", help = "Reference PXRD data file (CSV)")
parse.add_argument("ICDD_file", help = "ICDD data file (TXT)")

# ========== ADVANCED USER SETTINGS ==========
# Additional parameters for peak detection.
parse.add_argument("--intensity_threshold_value", type = float, default = 0.5,
                    help = "Intensity threshold for peak detection (default = 0.5)")
parse.add_argument("--angle_tolerance", type = float, default = 0.5,
                    help = "Angle tolerance for peak matching (default = 0.5)")
parse.add_argument("--intensity_tolerance", type = float, default = 0.1,
                    help = "Intensity tolerance for background filtering (default = 0.1)")    

args = parse.parse_args()
ref_file_name = args.ref_file
exp_file_name = args.exp_file
ICDD_file_name = args.ICDD_file
intensity_threshold_value = args.intensity_threshold_value
angle_tolerance = args.angle_tolerance
intensity_tolerance = args.intensity_tolerance 

# Physical Constants
"""
The chemical symbol of the incident radiation source (Cu or Mo) must be present in the experimental input file. 
The search is not case sensitive. 
"""
Wavelengths = {
    "Cu_Ka": 1.5406,  # Angstroms, Cu Kα radiation
    "Mo_Ka": 0.7093  # Angstroms, Mo Kα radiation    
    }

# Function for intensity normalisation: 
def intensity_norm(intensity):
    """
    Normalises intensity values relative to the maximum intensity. 

    Parameters
    ----------
    intensity : Python list of floats
        List of intensity values extracted from relevant data file

    Returns
    -------
    norm_intensity : Numpy array of floats
        Normalised intensity array
    """
    max_value = max(intensity)
    norm_intensity = np.asarray(intensity)/max_value
    return norm_intensity

# Reference data file key values

ref_two_theta = []
ref_intensity = []

# Parsing of the reference file 
with open(ref_file_name, 'r', encoding='latin-1') as ref_file: 
    next(ref_file)  # skip header line containing unnecessary information
    line = ref_file.readline()  
    while line:  # continue until file ends
    # turn the string into a list of strings separated at "," for csv file format
        split_line = line.split(",")
        ref_two_theta.append(float(split_line[0])) 
        ref_intensity.append(float(split_line[1])) 
        line = ref_file.readline()  # read next line
ref_normintensity = intensity_norm(ref_intensity)  

# Experimental data file key values
exp_two_theta = []
exp_intensity = []

# Parsing of the experimental data
with open(exp_file_name, 'r', encoding='latin-1') as exp_file:
    next(exp_file)  
    line = exp_file.readline()  
    while line:  
        split_line = line.split(",")
        exp_two_theta.append(float(split_line[0]))
        exp_intensity.append(float(split_line[1]))
        line = exp_file.readline()  
# Error check for corrupted experimental file
if not exp_two_theta or not exp_intensity:
    raise ValueError("Experimental data file appears empty or incorrectly formatted")
exp_normintensity = intensity_norm(exp_intensity)

# Finding the incident radiation wavelength 
with open(exp_file_name, 'r', encoding='latin-1') as exp_file:
    for line in exp_file:
        line_upper = line.upper() # makes sure it is not case sensitive
        if "CU" in line_upper:
            wavelength = Wavelengths["Cu_Ka"]
            break
        elif "MO" in line_upper:
            wavelength = Wavelengths["Mo_Ka"]
            break
    else:
        raise ValueError("Wavelength (Cu or Mo) not found in file")

# Title for the output section
print("-" * 58) # used for output layout purposes
print("THE FOLLOWING PROVIDES COMPREHENSIVE PEAK PICKING ANALYSIS")
print("-" * 58)
print("") # allows a space between various parts of the output

# Peak picking the most intense peaks
d_spacing_check = ppf.peak_pick(exp_two_theta, exp_normintensity, wavelength, intensity_threshold_value)
print("Most intense peaks:")
print("2θ value (°)  |  d spacing")
print("-" * 26)
for peak_positions, d_spacings in d_spacing_check:
    print(f"{peak_positions:12.2f}  |  {d_spacings:.2f}") # Values given to 2dp 
print("")

# Impurity peak picking
impurity_picker = ppf.impurity_pick(exp_two_theta, ref_two_theta, exp_normintensity, angle_tolerance, intensity_tolerance)
if impurity_picker: # only prints if impurities present
    print(f"There are impurity peaks at 2θ values of {impurity_picker}")
    print("")

# ICDD miller indices and associated diffraction angle
ICDD_two_theta = []
ICDD_miller_indices = []

# Parsing the ICDD file
with open(ICDD_file_name, 'r') as ICDD_file:
    next(ICDD_file)  
    for line in ICDD_file:
    # turns the string into a list of strings separated at "space" for txt file
        split_line = line.split() 
        if not split_line:  # Ignores lines with no/missing information
            continue
        ICDD_d_spacing = float(split_line[0])
        ICDD_theta_convert = math.asin(wavelength/(2*ICDD_d_spacing)) # using Bragg's Law
        ICDD_two_theta_convert = math.degrees(ICDD_theta_convert)*2
        ICDD_two_theta.append(ICDD_two_theta_convert)
        h = int(split_line[2])
        k = int(split_line[3])
        l = int(split_line[4].strip(",")) # An error occurs without stripping the , that appears attached to some of the l values
        ICDD_miller_indices.append((h,k,l))
# Defined zip for the 2θ values and miller indices from ICDD file, for printing the output
ICDD_peaks = list(zip(ICDD_two_theta, ICDD_miller_indices))
print("Complete peak annotation:")
print("2θ value (°)  |  Miller indices")
print("-" * 31)
for ICDD_angle, ICDD_miller in ICDD_peaks:
    if min(exp_two_theta) <= ICDD_angle <= max(exp_two_theta):
        print(f"{ICDD_angle:12.2f}  |  {ICDD_miller}")
print("")

# Plotting the overlayed PXRD pattern 
print("Plots:")
fig, ax = plt.subplots()

# Experimental pattern plot
ax.plot(exp_two_theta,exp_normintensity, linewidth = 0.5, label = "Experimental pattern", color = "black")
# Reference pattern plot
ax.plot(ref_two_theta,ref_normintensity, linewidth = 0.5, label = "Reference pattern", color = "blue")
# Ensuring a publication-quality plot
ax.set_title('Comparison PXRD Plot', fontsize = 14, fontweight = "bold")
plt.xlabel('2θ / °', fontsize = 12, fontweight = "bold")
plt.ylabel('Normalised intensity / a.u.', fontsize = 12, fontweight = "bold")
plt.legend()
# The x-axis is different for both plots of data so the axis is set to the experimental angle range
ax.set_xlim(min(exp_two_theta), max(exp_two_theta)) 
ax.set_ylim(0,1) # Intensities have been normalised from 0-1
# Save the plot
plt.savefig("Comparison_PXRD_Plot.png", dpi=300)
# Show the figure in the output
plt.show()

# Calling the Residual Analysis function
rmsd_analysis = ppf.rmsd(exp_two_theta, ref_two_theta, exp_normintensity, ref_normintensity, intensity_tolerance, angle_tolerance)
print(f"RMSD value is {rmsd_analysis[0]}")
print("")
print("Check for successful peak match:") 
print("Experimental 2θ value (°)  |  Reference 2θ value (°)")
print("-" * 52)
for experipeak, refpeak in rmsd_analysis[1]: # looks at the second return of the function
    print(f"{experipeak:25.2f}  |  {refpeak:.2f}")
    
