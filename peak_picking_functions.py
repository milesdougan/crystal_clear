# -*- coding: utf-8 -*-
"""
Created on Wed Feb 11 13:16:56 2026
@author: miles

Peak picking related functions, designed for usage within Crystal Clear however can be adapted for further programs.
 
"""
# Modules
import math
import numpy as np
import matplotlib.pyplot as plt


def peak_pick(two_theta, normalised_intensity, wavelength, intensity_threshold_value):
    """
    Identifies diffraction peaks and calculates d spacings for the most intense peaks using Bragg's Law.
    The function detects the local maxima in the normalised intensity data which exceeds an intensity threshold value. 

    Parameters
    ----------
    two_theta : Python list of floats
        The 2θ values experimentally recorded during PXRD
        Units: degrees (°)
    normalised_intensity : Numpy array of floats
        The intensity values normalised to 1, for the intensities experimentally recorded
    wavelength : Float from a dictionary
        The incident radiation wavelength value used in the experimental data recording
        Units: Angstroms (Å)
    intensity_threshold_value : Float
        The minimum intensity required for a peak to be picked

    Returns
    -------
    peak_check : a list of tuples or None
        A list of tuples displaying the diffraction angle and associated d spacing value, above a certain intensity.
        It is in the format 2θ value, d spacing.
        Units: degrees (°), Angstroms (Å)

    """
    peak_positions = []
    d_spacings = []
    for i in range(1, len(normalised_intensity) - 1): # Excludes endpoint values so the local maxima can be detected using neighbouring values.
        if (normalised_intensity[i] >= intensity_threshold_value
            and normalised_intensity[i] > normalised_intensity[i + 1]
            and normalised_intensity [i] > normalised_intensity[i - 1]):
            theta_radians = math.radians(two_theta[i]/2)
            d_spacing = wavelength / (2*math.sin(theta_radians))
            d_spacings.append(d_spacing)
            peak_positions.append(two_theta[i])
    if peak_positions: # make sure the loop iterates through all peaks.
        peak_check = list(zip(peak_positions, d_spacings))    
        return peak_check
    else:
        return None
    

def impurity_pick(exp_two_theta, ref_two_theta, exp_normintensity, angle_tolerance, intensity_tolerance):    
    """
    Picks out the impurity peaks in an experimental PXRD dataset.
    The function follows a similar convention to that of peak_pick:
    It detects the local maxima in the normalised intensity experimental data which exceeds an intensity cutoff.
    The detected experimental peaks are then compared to reference peaks to elucidate impurities outside an angle tolerance.

    Parameters
    ----------
    exp_two_theta : Python list of floats
        The 2θ values/diffraction angles experimentally recorded during PXRD
        Units: degrees (°)
    ref_two_theta : Python list of floats
        The 2θ values/diffraction angles for the peaks in the reference data
        Units: degrees (°)
    exp_normintensity : NumPy array of floats
        The intensity values normalised to 1, for the intensities experimentally recorded
    angle_tolerance : float
        The allowed diffraction angle difference between experimental and reference peaks for them to be a match
        Units: degrees (°)
    intensity_tolerance : float
        The cutoff intensity for distinguishing a peak from the background intensity.

    Returns
    -------
    impurity_peaks : python list of floats
        A list of 2θ values that correspond to impurity peaks present in the experimental data.
        Units: degrees (°)
        
    """
    impurity_peaks = [] # 2theta values
    for i in range(1, len(exp_normintensity) - 1): # Excludes endpoint values so the local maxima can be detected using neighbouring values. 
        if (exp_normintensity[i] >= intensity_tolerance
            and exp_normintensity[i] > exp_normintensity[i+1] 
            and exp_normintensity[i] > exp_normintensity[i-1]):
            exp_peak = exp_two_theta[i]
            peak_found = False
            for ref_peak in ref_two_theta:
                if abs(exp_peak - ref_peak) <= angle_tolerance: # peak present in reference
                       peak_found = True
                       break
            if peak_found == False:
                impurity_peaks.append(exp_peak)
    return impurity_peaks


def rmsd(exp_two_theta, ref_two_theta, exp_normintensity, ref_normintensity, intensity_tolerance, angle_tolerance):
    """
    Calculates the root mean squared deviation (RMSD) between experimental and reference 2θ peak positions.  
    It detects the local maxima in the experimental and reference normalised intensity data which exceeds an intensity cutoff.
    The detected experimental and reference peaks are then compared to identify peak match within an angle tolerance.
    This angular difference contributes to the RMSD calculation.
    If a detected peak is not matched, an impurity peak is highlighted, also contributing to the RMSD value. 
    Background intensity is disregarding.

    Parameters
    ----------
    exp_two_theta, ref_two_theta : Python list of floats
        Diffraction angles/ 2θ value 
        Units: degrees (°)
    exp_normintensity, ref_normintensity : Numpy array of floats
        The intensity values normalised to 1
    intensity_tolerance : float
        The cutoff intensity for distinguishing a peak from the background intensity.
    angle_tolerance : float
        The allowed diffraction angle difference between experimental and reference peaks for them to be a match.
        Units: degrees (°)
    
    Returns
    -------
    rmsd_out : float
        A value for the Root Mean Square Deviation (RMSD) between the experimental and reference peak positions.
    compare_peaks : list of tuples
        A list of the matched peaks in the format Experimental 2θ value, Reference 2θ value 
        Units: degrees (°), degrees (°)
    plt : plot
        A plot of peak position differences (Δ2θ) versus experimental 2θ.
    ValueError : raise
        Causes an error if no peaks are matched between the experimental and reference data.
        
    """
    experipeak = []
    extract_impurity = []
    refpeak = []
    ref_angles = []
    sum_squares = 0
    peaks_matched = 0
    
    # Find reference peak maxima
    for i in range(1, len(ref_normintensity) - 1):
        if (ref_normintensity[i] >= intensity_tolerance
            and ref_normintensity[i] > ref_normintensity[i+1] 
            and ref_normintensity[i] > ref_normintensity[i-1]):
            ref_angles.append(ref_two_theta[i])
    ref_angles = np.array(ref_angles)
    
    # Find experimental peak maxima
    for i in range(1, len(exp_normintensity) - 1):
        if (exp_normintensity[i] >= intensity_tolerance
            and exp_normintensity[i] > exp_normintensity[i+1] 
            and exp_normintensity[i] > exp_normintensity[i-1]):
            exp_angle = exp_two_theta[i]
            
            # Difference in peak matching
            peak_dif = np.abs(exp_angle - ref_angles)
            if len(peak_dif) == 0: # if no reference peak is found
                extract_impurity.append(round(exp_angle, 2))
                sum_squares += exp_angle**2 # impurity added to RMSD
                continue
            # Find the closest peak match
            min_index = np.argmin(peak_dif)             
            if peak_dif[min_index] <= angle_tolerance:
                ref_angle = ref_angles[min_index]
                experipeak.append(round(exp_angle, 2))
                refpeak.append(round(ref_angle, 2))
                sum_squares += (exp_angle - ref_angle)**2 # peak match added to RMSD
                peaks_matched += 1 # For normalisation against number of peaks matched
            else:
                extract_impurity.append(round(exp_angle, 2))
                sum_squares += exp_angle**2
    # if no peaks match the RMSD cannot be calculated
    if peaks_matched == 0 :
        raise ValueError("RMSD cannot be computed due to zero matched peaks")
                
    rmsd_out = np.sqrt(sum_squares/peaks_matched)
    rmsd_out = np.round(rmsd_out, 4) 
    compare_peaks = list(zip(experipeak, refpeak)) # matched peak pairs checked for mismatch
    
    # RMSD Plot
    if rmsd_out:
        experipeak = np.array(experipeak)
        refpeak = np.array(refpeak)
        delta_theta = abs(experipeak - refpeak)
        
    fig, ax = plt.subplots()
    # Plot of Δ2θ Vs. 2θ
    ax.plot(experipeak, delta_theta, "x--", linewidth = 0.5, color = "r")
    plt.xlabel('2θ / °', fontsize = 12, fontweight = "bold" )
    plt.ylabel('Δ2θ / °', fontsize = 12, fontweight = "bold")
    ax.set_title('Peak Differences Vs. 2θ', fontsize = 14, fontweight = "bold")    
    # Save plot
    plt.savefig("RMSD_Peak_Difference_Plot.png", dpi=300)
    # Displays the figure in the output
    plt.show()
    
    return rmsd_out, compare_peaks

