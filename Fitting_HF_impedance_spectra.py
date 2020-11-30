#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Fitting_HF_impedance_spectra.py, version 0.1.

This script is for fitting high-frequency linear impedance branches of
PEMFC and PEMWE impedance spectra with analytical models from references [1,2].
Effective electronic and protonic conductivities are extracted, and
corresponding effective ohmic resistances are computed.

Please refer to README.rst for documentation.

Requirements: Python 3 with the following packages:
              matplotlib, numpy, scipy, sympy, time, math, csv, cmath, mpmath,
              sys, getopt

References
[1] A. Kulikovsky. "A model for impedance of a PEM fuel cell cathode with
    poor electron conductivity." Journal of Electroanalytical Chemistry 801
    (2017): 122-128.
[2] A. Kosakian, M. Secanell. "Estimating charge-transport properties of 
    fuel-cell and electrolyzer catalyst layers via electrochemical impedance 
    spectroscopy." Submitted to Electrochimica Acta (2020).

Copyright 2020 Aslan Kosakian, Marc Secanell

               Energy Systems Design Laboratory
               University of Alberta
               Canada

Distributed under the MIT license. See LICENSE.txt for details.
"""

#################
# Initial setup #
#################

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
from scipy import stats
from scipy import optimize
from sympy import solve
from sympy import Symbol
from matplotlib import colors
import time
import math
import csv
import cmath
import mpmath as mp
import sys, getopt

############################
# Reading and writing data #
############################

# Reads the parameter file. Ignores comments (start with #).
# Parameter format: "parameter = value".
def read_parameter_file(filename):
    contents =  open(filename, 'r')
    lines = contents.readlines()
    dict_data = {}
    for line in lines:
        strip_split_line = line.strip().split('#',1)[0]
        if strip_split_line!='':
            prm_line_split = strip_split_line.split('=')
            dict_data.update( {prm_line_split[0].strip() : prm_line_split[1]} )
    return dict_data
        
        
# Reads frequency and impedance data from a text file (any text format).
# The format of the file must be as follows.
# The first line contains headers that are named "Frequency_Hz",
# "Z_re_mOhm*cm2", and "Z_neg_imag_mOhm*cm2" (no quotation marks).
# The second and subsequent rows contain data corresponding to the
# headers. The delimeter is determined automatically.
def load_impedance_from_file(filename):
    with open(filename, newline='') as csvfile:
        # Automatically determine the delimiter
        dialect = csv.Sniffer().sniff(csvfile.read(1024))
        # Return to the beginning of the file
        csvfile.seek(0)
        # Read the contents of the file
        contents = csv.reader(csvfile, dialect)
        # Process the file row by row
        frequency_data = []
        Z_re_data = []
        Z_neg_imag_data = []
        for row in contents:
            if contents.line_num>1:
                frequency_data.append(float(row[0]))
                Z_re_data.append(float(row[1]))
                Z_neg_imag_data.append(float(row[2]))
    return(zip(*sorted(zip(frequency_data,Z_re_data,Z_neg_imag_data))))

# Writes requency and impedance data to a text file.
# The format of the file is as follows.
# The first line contains headers that are named "Frequency_Hz",
# "Z_re_mOhm*cm2", and "Z_neg_imag_mOhm*cm2" (no quotation marks).
# The second and subsequent rows contain data corresponding to the
# headers. Tab symbol ("\t") is used as a delimiter.
def save_impedance_to_file(filename,frequencies,Z_real,Z_neg_imag):
    output = open(filename,"w")
    output.write("Frequency_Hz" + "\t" + "Z_re_mOhm*cm2" + "\t" + "Z_neg_imag_mOhm*cm2" + "\n")
    for i in range(0,len(frequencies)):
        output.write(str(frequencies[i]) + "\t" + str(Z_real[i]) + "\t" + str(Z_neg_imag[i]) + "\n")

#####################
# Analytical models #
#####################

# This is equation (36) from reference [1] in dimensional form (for H2/O2 
# impedance spectra and for H2/N2 spectra with significant hydrogen 
# cross-over; see discussion in reference [2]).
#
# The original variable notation from reference [1] is used.
#
# Input arguments and units:
#   effective electronic conductivity (sigma_e), S/cm;
#   effective protonic conductivity (sigma_p), S/cm;
#    catalyst-layer thickness (L_CL) in cm;
#   volumetric double-layer capacitance (Cdl_vol), F/cm3;
#   signal frequency (frequency), Hz;
#   operating current density (j0), A/cm2;
#   Tafel slope (b), V.
#
# Returns analytical impedance (a complex number).
def analytical_H2_O2_impedance_Kulikovsky(sigma_e,sigma_p,L_CL,Cdl_vol,frequency,j0,b):
    omega = 2*np.pi*frequency
    k_sigma = sigma_e/(sigma_p+1e-14)
    Omega_tilde = omega*Cdl_vol*L_CL**2/(sigma_p+1e-14)
    j0_tilde = j0*L_CL/(sigma_p*b+1e-14)
    p0 = cmath.sqrt(-(j0_tilde+1j*Omega_tilde)*(1+1.0/k_sigma))
    q0 = cmath.sqrt(2*k_sigma*(1+k_sigma)*(1j*j0_tilde-Omega_tilde))
    tmp1 = complex(mp.fmul(1j*q0,complex(mp.sin(p0))))
    tmp2 = complex(mp.fmul((1+k_sigma**2),complex(mp.cos(p0))))
    tmp3 = complex(mp.fmul((1j*q0-(1+1j)*p0),complex(mp.sin(p0))))
    tmp31 = complex(mp.fadd(2*k_sigma,tmp2))
    tmp32 = complex(mp.fmul(1+1j,tmp31))
    tmp33 = complex(mp.fmul(k_sigma,tmp3))
    tmp4 = complex(mp.fdiv((tmp1+tmp32),(tmp33+1e-14)))
    Z_Kulikovsky = 1000*L_CL/sigma_p*tmp4
    return Z_Kulikovsky

# Vectorized version of analytical_H2_O2_impedance_Kulikovsky
v_analytical_H2_O2_impedance_Kulikovsky = np.vectorize(analytical_H2_O2_impedance_Kulikovsky)

# Returns the residual between the input data and the analytical impedance
# computed with analytical_H2_O2_impedance_Kulikovsky. Three residual methods
# are supported: absolute difference (magnitude of Z_input-Z_model, "abs"),
# relative difference (absolute difference normalized by Z_model, "rel"), and
# absolute difference normalized by the magnitude of Z_model ("rel_mag"). The
# final residual value is divided by the square root of its length (in analogy
# with root-mean-square deviation).
#
# Input parameters (other than those for analytical_H2_O2_impedance_Kulikovsky)
# are:
#   a tuple of electronic and protonic conductivities (x), all in S/cm;
#   input real impedance (Z_re_input), mOhm*cm2;
#   input negative imaginary impedance (Z_neg_imag_input), mOhm*cm2;
#   residual method (method);
#   a verbose parameter to print the residual (verbose).
#
# Note that a vector of frequencies is supplied in this case (Hz).
def residual_analytical_H2_O2_impedance_Kulikovsky(x,Z_re_input,Z_neg_imag_input,L_CL,Cdl_vol,frequencies,j0,b,method,verbose=False):
    sigma_e = x[0]
    sigma_p = x[1]
    residual = 0
    
    for i in range(0,len(frequencies)):
        Z_Kulikovsky = analytical_H2_O2_impedance_Kulikovsky(sigma_e,sigma_p,L_CL,Cdl_vol,frequencies[i],j0,b)
        if method=="abs":
            residual = residual + ((Z_re_input[i]-np.real(Z_Kulikovsky)))**2 + ((-Z_neg_imag_input[i]-np.imag(Z_Kulikovsky)))**2
        elif method=="rel":
            residual = residual + ((Z_re_input[i]-np.real(Z_Kulikovsky))/(np.real(Z_Kulikovsky)+1e-14))**2 + ((-Z_neg_imag_input[i]-np.imag(Z_Kulikovsky))/(np.imag(Z_Kulikovsky)+1e-14))**2
        elif method=="rel_mag":
            Z_mag = np.sqrt(Z_re_input[i]**2+Z_neg_imag_input[i]**2)
            residual = residual + ((Z_re_input[i]-np.real(Z_Kulikovsky))/(np.abs(Z_Kulikovsky)+1e-14))**2 + ((-Z_neg_imag_input[i]-np.imag(Z_Kulikovsky))/(np.abs(Z_Kulikovsky)+1e-14))**2
        else:
            print("[Error] Unsupported residual method.")
            sys.exit(1)
        
    residual = np.sqrt(residual/len(Z_re_input))
    if math.isnan(residual):
        # If the analytical model returns NaN (for example, due to an
        # overflow in some trigonometric functions), replace it with
        # a large number for the fitting algorithm to work. This is not
        # optimal for gradient-based optimizationalgorithms, but they are
        # not used in this code.
        residual = float("inf")
    if verbose:
        print("residual="+str(residual)+" for sigma_e="+str(sigma_e)+" S/cm, sigma_p="+str(sigma_p)+" S/cm.")
    return residual


# This is equation (44) from reference [1] in dimensional form (for H2/N2 
# spectra without hydrogen cross-over; see discussion in reference [2]).
#
# The original variable notation from reference [1] is used.
#
# Input arguments and units:
#   effective electronic conductivity (sigma_e), S/cm;
#   effective protonic conductivity (sigma_p), S/cm;
#   catalyst-layer thickness (L_CL) in cm;
#   volumetric double-layer capacitance (Cdl_vol), F/cm3;
#   signal frequency (frequency), Hz.
#
# Returns analytical impedance (a complex number)
def analytical_H2_N2_impedance_Kulikovsky(sigma_e,sigma_p,L_CL,Cdl_vol,frequency):
    omega = 2*np.pi*frequency
    k_sigma = sigma_e/(sigma_p+1e-14)
    Omega_tilde = omega*Cdl_vol*L_CL**2/(sigma_p+1e-14)
    p0 = cmath.sqrt(-1j*Omega_tilde*(1+1.0/k_sigma))
    q0 = cmath.sqrt(-2*k_sigma*(1+k_sigma)*Omega_tilde)
    tmp1 = complex(mp.fmul(1j*q0,complex(mp.sin(p0))))
    tmp2 = complex(mp.fmul((1+k_sigma**2),complex(mp.cos(p0))))
    tmp3 = complex(mp.fmul((1j*q0-(1+1j)*p0),complex(mp.sin(p0))))
    tmp31 = complex(mp.fadd(2*k_sigma,tmp2))
    tmp32 = complex(mp.fmul(1+1j,tmp31))
    tmp33 = complex(mp.fmul(k_sigma,tmp3))
    tmp4 = complex(mp.fdiv((tmp1+tmp32),(tmp33+1e-14)))
    Z_Kulikovsky = 1000*L_CL/sigma_p*tmp4
    return Z_Kulikovsky

# Vectorized version of analytical_H2_N2_impedance_Kulikovsky
v_analytical_H2_N2_impedance_Kulikovsky = np.vectorize(analytical_H2_N2_impedance_Kulikovsky)

# Returns the residual between the input data and the analytical impedance
# computed with analytical_H2_N2_impedance_Kulikovsky. Three residual methods
# are supported: absolute difference (magnitude of Z_input-Z_model, "abs"),
# relative difference (absolute difference normalized by Z_model, "rel"), and
# absolute difference normalized by the magnitude of Z_model ("rel_mag"). The
# final residual value is divided by the square root of its length (in analogy
# with root-mean-square deviation).
#
# Input parameters (other than those for analytical_H2_N2_impedance_Kulikovsky)
# are:
#   a tuple of electronic and protonic conductivities (x), all in S/cm;
#   input real impedance (Z_re_input), mOhm*cm2;
#   input negative imaginary impedance (Z_neg_imag_input), mOhm*cm2;
#   residual method (method);
#   a verbose parameter to print the residual (verbose).
#
# Note that a vector of frequencies is supplied in this case (Hz).
def residual_analytical_H2_N2_impedance_Kulikovsky(x,Z_re_input,Z_neg_imag_input,L_CL,Cdl_vol,frequencies,method,verbose=False):
    sigma_e = x[0]
    sigma_p = x[1]
    residual = 0
    
    for i in range(0,len(frequencies)):
        Z_Kulikovsky = analytical_H2_N2_impedance_Kulikovsky(sigma_e,sigma_p,L_CL,Cdl_vol,frequencies[i])
        if method=="abs":
            residual = residual + ((Z_re_input[i]-np.real(Z_Kulikovsky)))**2 + ((-Z_neg_imag_input[i]-np.imag(Z_Kulikovsky)))**2
        elif method=="rel":
            residual = residual + ((Z_re_input[i]-np.real(Z_Kulikovsky))/(np.real(Z_Kulikovsky)+1e-14))**2 + ((-Z_neg_imag_input[i]-np.imag(Z_Kulikovsky))/(np.imag(Z_Kulikovsky)+1e-14))**2
        elif method=="rel_mag":
            Z_mag = np.sqrt(Z_re_input[i]**2+Z_neg_imag_input[i]**2)
            residual = residual + ((Z_re_input[i]-np.real(Z_Kulikovsky))/(np.abs(Z_Kulikovsky)+1e-14))**2 + ((-Z_neg_imag_input[i]-np.imag(Z_Kulikovsky))/(np.abs(Z_Kulikovsky)+1e-14))**2
        else:
            print("[Error] Unsupported residual method.")
            sys.exit(1)
    
    residual = np.sqrt(residual/len(Z_re_input))
    if math.isnan(residual):
        # If the analytical model returns NaN (for example, due to an
        # overflow in some trigonometric functions), replace it with
        # a large number for the fitting algorithm to work. This is not
        # optimal for gradient-based optimizationalgorithms, but they are
        # not used in this code.
        residual = float("inf")
    if verbose:
        print("residual="+str(residual)+" for s_e="+str(sigma_e)+" S/cm, s_p="+str(sigma_p)+" S/cm.")
    return residual

# Returns effective ohmic resistance as R=L_CL/(3 sigma), where L_CL is
# the CL thickness and sigma is an effective conductivity
def compute_ohmic_resistance(L_CL,sigma):
    return L_CL/3.0/sigma*1000

##############################
# Data processing and output #
##############################

# Identifies the linear high-frequency impedance branch for fitting. This is
# done by comparing the HFR-corrected real impedance and negative imaginary
# impedance with a given relative tolerance.
#
# Input parameters are:
#   frequency (f), Hz;
#   input real impedance (Z_re_input), mOhm*cm2;
#   input negative imaginary impedance (Z_neg_imag_input), mOhm*cm2;
#   approximate value of Z_re_input where the linear branch ends (fit_max_Z_re_input), mOhm*cm2;
#   relative threshold that will be used to dentify the linear branch (threshold).
#
# Returns the truncated impedance (real and negative imaginary parts) and
# frequency data.
def select_HF_linear_impedance_branch(f,Z_re_input,Z_neg_imag_input,fit_max_Z_re_input,threshold):
    # Form an HFR-corrected real impedance (necessary for the correct
    # determination of the 45-degree branch)
    Z_re_input_corr = [x-min(Z_re_input) for x in Z_re_input]
    index_range_45_branch = []
    for i in range(0,len(Z_re_input_corr)):
       if abs(Z_re_input_corr[i]-Z_neg_imag_input[i])/abs(Z_re_input_corr[i]+1e-14)<threshold and abs(Z_re_input[i]+1e-14)<fit_max_Z_re_input:
           index_range_45_branch.append(i)
    if len(index_range_45_branch)==0:
        print("[Error] High-frequency spectrum is not linear within the given tolerance. No data points were identified.")
        sys.exit(1)
    if len(index_range_45_branch)<10:
        print("[Warning] Only "+str(len(index_range_45_branch))+ \
              " were identified as a part of the linear impedance branch at high frequencies \
                within the given tolerance. The fitted results may be inaccurate.")
    return Z_re_input[min(index_range_45_branch):max(index_range_45_branch)],  \
           Z_neg_imag_input[min(index_range_45_branch):max(index_range_45_branch)], \
           f[min(index_range_45_branch):max(index_range_45_branch)]
       
# Performs fitting with the H2/O2 and H2/N2 methods (non-graphical).
#
# Input parameters are:
#   residual function to be called (res_func);
#   a tuple of conductivity ranges for computing the residual (s_ranges);
#   a tuple of arguments for the residal function (args);
#   number of CPU threads to be used for fitting (n_t);
#   input real impedance (Z_re_in), mOhm*cm2;
#   input negative imaginary impedance (Z_neg_imag_in), mOhm*cm2.
#
# Returns:
#   final residual (res);
#   fitted conductivities (s_e, s_p), S/cm;
#   fitted real and negative imaginary impedance (Z_re_fit, Z_neg_imag_fit), mOhm*cm2;
#   R^2 values of the fits for the real and imaginary parts (R2_re, R2_neg_imag).
def perform_curve_fitting(res_func,s_ranges,args,n_t,fitting_model):
    Z_re_in = args[0]
    Z_neg_imag_in = args[1]
    L_CL = args[2]
    Cdl_cm3 = args[3]
    f_in = args[4]
    residual_method = args[-2]
    
    start = time.time()
    res = optimize.brute(res_func, \
                         ranges=s_ranges,\
                         args=args,\
                         workers=n_t,\
                         full_output=True)
    print("Fitting finished in {:.1f} seconds.".format((time.time() - start)))
    str_units = ""
    if residual_method=="abs":
        str_units = " mOhm*cm2"
    print("Residual: "+str(res[1])+str_units)
    print("")
        
    s_e = res[0][0]
    s_p = res[0][1]
    
    if fitting_model=="H2_O2":
        dc_current_density = args[-4]
        Tafel_slope = args[-3]
        Z_fit = v_analytical_H2_O2_impedance_Kulikovsky(s_e,s_p,L_CL,Cdl_cm3,f_in,dc_current_density,Tafel_slope)
    elif fitting_model=="H2_N2":
        Z_fit = v_analytical_H2_N2_impedance_Kulikovsky(s_e,s_p,L_CL,Cdl_cm3,f_in)
    else:
        print("[Error] Unsupported fitting method passed to perform_curve_fitting.")
        sys.exit(1)
    Z_re_fit = np.real(Z_fit)
    Z_neg_imag_fit = [-x for x in np.imag(Z_fit)]
    
    slope, intercept, R2_re, p_value, std_err = stats.linregress(Z_re_in,Z_re_fit)
    slope, intercept, R2_neg_imag, p_value, std_err = stats.linregress(Z_neg_imag_in,Z_neg_imag_fit)
        
    return res,s_e,s_p,Z_re_fit,Z_neg_imag_fit,R2_re,R2_neg_imag

# Outputs fitted conductivities, computed resistances, and R^2 values.
#
# Inputs:
#   fitted conductivities (sigma_e, sigma_p), S/cm;
#   R^2 values (R2_1, R2_2), the meaning of which depends on the fitting method;
#   CL thickness (L_CL), cm;
#   fitting model ("H2_O2","H2_N2","H2_N2_graphical").
def output_fitting_info(sigma_e,sigma_p,R2_1,R2_2,L_CL,fitting_model):
    print("Fitted effective electronic conductivity: ",sigma_e," S/cm.")
    print("Fitted effective protonic conductivity: ",sigma_p," S/cm.")
    print("Corresponding effective electronic resistance: ",compute_ohmic_resistance(L_CL,sigma_e)," mOhm*cm2.")
    print("Corresponding effective protonic resistance: ",compute_ohmic_resistance(L_CL,sigma_p)," mOhm*cm2.")
    if fitting_model=="H2_O2" or fitting_model=="H2_N2":
        print("R^2: ",R2_1," (real part), ",R2_2," (negative imaginary part).")
    elif fitting_model=="H2_N2_graphical":
        print("R^2: ",R2_1," (linear low-frequency fit), ",R2_2," (linear high-freqiency fit).")
    print("")
    msg = "[Note] If conductivity (resistance) A is larger than conductivity (resistance) B against " + \
          "the physical expectation,  swap the results. Swapped conductivities (resistances) are equally valid."
    print(msg)

# Checks whether the input is a positive float.
def process_input_float(val,name):
    try:
        if float(val)<-1e-14:
            print("[Input error] Value of "+name+" must be non-negative ('"+name+"=').")
            sys.exit(1)
    except ValueError:
        print("[Input error] Value of "+name+" must be a float.")
        sys.exit(1)
    return float(val)

# Checks whether the input is a positive integer.
def process_input_int(val,name):
    try:
        if int(val)<1:
            print("[Input error] Value of "+name+" must be positive ('"+name+"=').")
            sys.exit(1)
    except ValueError:
        print("[Input error] Value of "+name+" must be an integer.")
        sys.exit(1)
    return int(val)

# Checks whether the input is a bolean.
def process_input_bool(val,name):
    if val=="True":
        return True
    elif val=="False":
        return False
    else:
        print("[Input error] Value of "+name+" must be a boolean.")
        sys.exit(1)

# Check if the requsted parameter was read from the parameter file.
# Inputs are dict_prm, a parameter dictionary read from the input parameter
# file by read_parameter_file function, and a parameter key.
def parse_parameter(dict_prm,key,prm_type):
    try:
        if prm_type=="float":
            value = process_input_float(dict_prm[key].strip(),key)
        elif prm_type=="int":
            value = process_input_int(dict_prm[key].strip(),key)
        elif prm_type=="bool":
            value = process_input_bool(dict_prm[key].strip(),key)
        elif prm_type=="string":
            value = dict_prm[key].strip()
        else:
            print("[Error] Unsupported prm_type passed to parse_parameter function.")
            sys.exit(1)
    except KeyError:
        print("[Error] Parameter '" + str(key) + "' not found in the parameter file.")
        print("        This error may occur when you either have a typo in the")
        print("        parameter name or forgot to include a parameter required")
        print("        for the chosen fitting method.")
        sys.exit(1)
    return value

########
# Main #
########

def main(argv):
    ##############
    # Data input #
    ##############
    
    prm_filename = ''
    verbose = False
    n_threads=1
    plot_input_data = False
    
    try:
        opts, args = getopt.getopt(argv,"hv",["help","prm=","threads=","verbose","plot_input_data"])
    except getopt.GetoptError:
        print('[Error] Unrecognized option(s).')
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            print("Fitting_HF_impedance_spectra.py, version 0.1")
            print("Copyright 2020 Aslan Kosakian, Marc Secanell")
            print("Energy Systems Design Laboratory, University of Alberta, Canada")
            print("")
            print("Distributed under the MIT license. See LICENSE.txt for details.")
            print("")
            print("Options:")
            print("-h, --help            Output this message")
            print("--prm=                Input parameter file")
            print("--threads=            Number of CPU threads to use")
            print("--verbose             Output fitting residual at each iteration")
            print("--plot_input_data     Plot the input data")
            print("")
            print("Please refer to README.rst for documentation.")
            sys.exit()
        elif opt in ("-v", "--verbose"):
            verbose = True
        elif opt=="--prm":
            prm_filename = arg
        elif opt=="--threads":
            n_threads = process_input_int(arg,"threads")
        elif opt=="--plot_input_data":
            plot_input_data = True
       
    if prm_filename=='':
        print("[Error] No parameter file provided. Use the '--prm=' flag.")
        sys.exit(1)
    
    dict_prm                         = read_parameter_file(prm_filename)
    filename                         = parse_parameter(dict_prm,"Input file","string")
    output_filename                  = parse_parameter(dict_prm,"Output file","string")
    fitting_model                    = parse_parameter(dict_prm,"Fitting method","string")
    Cdl_cm3             = parse_parameter(dict_prm,"Volumetric double-layer capacitance, F/cm3","float")
    L_CL                             = parse_parameter(dict_prm,"Catalyst-layer thickness, um","float")/1e4 # um to cm
    fit_max_re_Z                     = parse_parameter(dict_prm,"Real impedance cutoff, mOhm*cm2","float")
    threshold                        = parse_parameter(dict_prm,"Linear branch threshold","float")
    HFR_correction                   = parse_parameter(dict_prm,"HFR correction, mOhm*cm2","float")
    Bode_log_scale                   = parse_parameter(dict_prm,"Use log scale in Bode plot","bool")
    if fitting_model in ("H2_O2","H2_N2"):
        sigma_e_search_min  = parse_parameter(dict_prm,"Minimum electronic conductivity, S/cm","float")
        sigma_e_search_max  = parse_parameter(dict_prm,"Maximum electronic conductivity, S/cm","float")
        sigma_e_search_step = parse_parameter(dict_prm,"Electronic-conductivity step size, S/cm","float")
        sigma_p_search_min  = parse_parameter(dict_prm,"Minimum protonic conductivity, S/cm","float")
        sigma_p_search_max  = parse_parameter(dict_prm,"Maximum protonic conductivity, S/cm","float")
        sigma_p_search_step = parse_parameter(dict_prm,"Protonic-conductivity step size, S/cm","float")
        residual_method     = parse_parameter(dict_prm,"Residual method","string")
        output_full_analytical_impedance = parse_parameter(dict_prm,"Output whole fitted spectrum","bool")
        if fitting_model=="H2_O2":
            dc_current_density = parse_parameter(dict_prm,"Current density, A/cm2","float")
            Tafel_slope        = parse_parameter(dict_prm,"Tafel slope, V","float")
    elif fitting_model=="H2_N2_graphical":
        n_LF_points_to_use = parse_parameter(dict_prm,"Number of low-frequency points to fit","int")
    else:
        print("[Error] Unsupported fitting method '" + str(fitting_model) + "'.")
        sys.exit(1)
       
    print("Reading data from "+filename+"...")
    frequencies_input, Z_re_input, Z_neg_imag_input = load_impedance_from_file(filename) 

    #################
    # Plot settings #
    #################
    
    dpi_val=80.0
    Fontsize = 16
    plt.rc('font', size=Fontsize-3)          # controls default text sizes
    plt.rc('axes', titlesize=Fontsize)     # fontsize of the axes title
    plt.rc('axes', labelsize=Fontsize)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=Fontsize)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=Fontsize)    # fontsize of the tick labels
    plt.rc('legend', fontsize=Fontsize-2)    # legend fontsize

    # If the user just wants to plot the input data, do so and exit
    if plot_input_data:
        # ------------ #
        # Nyquist plot #
        # ------------ #
        fig1 = plt.figure(1,facecolor='white', figsize=(500/dpi_val, 400/dpi_val), dpi=dpi_val)
        
        ax1 = plt.subplot2grid((1, 1), (0, 0))
        
        ax1.plot(Z_re_input,Z_neg_imag_input,marker="o",fillstyle="none",color="k",linestyle="none",markersize=10,label='Input data', zorder=5)
        
        ax1.set_xlabel(r'Re(Z), mOhm*cm2')
        ax1.set_ylabel(r'-Im(Z), mOhm*cm2')
        
        ax1.set_aspect('equal')
        
        ax1.legend(loc='best')
        
        fig1.tight_layout()
        #plt.show()
        
        # --------- #
        # Bode plot #
        # --------- #
        
        fig2 = plt.figure(2,facecolor='white', figsize=(500/dpi_val, 400/dpi_val), dpi=dpi_val)
        
        ax2 = plt.subplot2grid((1, 1), (0, 0))
        
        ax2.plot(frequencies_input,Z_neg_imag_input,marker="o",fillstyle="none",color="k",linestyle="none",markersize=10,label='Original data', zorder=5)
        
        ax2.set_xlabel(r'Frequency_Hz')
        ax2.set_ylabel(r'-Im(Z), mOhm*cm2')
        
        if Bode_log_scale:
            ax2.set_xscale('log')
            locmaj = ticker.LogLocator(base=10,numticks=15) 
            ax2.xaxis.set_major_locator(locmaj)
            locmin = ticker.LogLocator(base=10.0,subs=(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),numticks=12)
            ax2.xaxis.set_minor_locator(locmin)
            ax2.xaxis.set_minor_formatter(ticker.NullFormatter())
        
        ax2.legend(loc='best')
        
        fig2.tight_layout(pad=1.5)
        plt.show()
        sys.exit(0)

    ################
    # Data fitting #
    ################

    # HFR correction
    Z_re_input = np.array(Z_re_input) - HFR_correction
    fit_max_re_Z = fit_max_re_Z - HFR_correction

    Z_re_input_trunc, Z_neg_imag_input_trunc, frequencies_input_trunc = select_HF_linear_impedance_branch(frequencies_input, Z_re_input, Z_neg_imag_input, fit_max_re_Z, threshold)
    print("Successfully identified the linear high-frequency impedance branch; " + str(len(frequencies_input_trunc)) + " data points will be used for fitting.")
    
    if fitting_model in ("H2_O2","H2_N2"):
        sigma_ranges = (slice(sigma_e_search_min,sigma_e_search_max,sigma_e_search_step), slice(sigma_p_search_min,sigma_p_search_max,sigma_p_search_step))
        arguments_base = (Z_re_input_trunc,Z_neg_imag_input_trunc,L_CL,Cdl_cm3,frequencies_input_trunc)
    
        if fitting_model=="H2_O2":
            print("Fitting the data with the H2/O2 model...")
            arguments_H2_O2 = arguments_base + (dc_current_density,Tafel_slope,residual_method,verbose,)
            residual,sigma_e_fit,sigma_p_fit,Z_re_fit,Z_neg_imag_fit,r_value_re,r_value_neg_imag = perform_curve_fitting(residual_analytical_H2_O2_impedance_Kulikovsky,
                                                                                                                         sigma_ranges,
                                                                                                                         arguments_H2_O2,
                                                                                                                         n_threads,
                                                                                                                         fitting_model)
            output_fitting_info(sigma_e_fit,sigma_p_fit,r_value_re,r_value_neg_imag,L_CL,fitting_model)
            if output_full_analytical_impedance:
                Z_fit_full = v_analytical_H2_O2_impedance_Kulikovsky(sigma_e_fit,sigma_p_fit,L_CL,Cdl_cm3,frequencies_input,dc_current_density,Tafel_slope)
                save_impedance_to_file(output_filename,frequencies_input,np.real(Z_fit_full),[-x for x in np.imag(Z_fit_full)])
            else:
                save_impedance_to_file(output_filename,frequencies_input_trunc,Z_re_fit,Z_neg_imag_fit)
        elif fitting_model=="H2_N2":
            print("Fitting the data with the H2/N2 model...")
            arguments_H2_N2 = arguments_base + (residual_method,verbose,)
            residual,sigma_e_fit,sigma_p_fit,Z_re_fit,Z_neg_imag_fit,r_value_re,r_value_neg_imag = perform_curve_fitting(residual_analytical_H2_N2_impedance_Kulikovsky,
                                                                                                                         sigma_ranges,
                                                                                                                         arguments_H2_N2,
                                                                                                                         n_threads,
                                                                                                                         fitting_model)
            output_fitting_info(sigma_e_fit,sigma_p_fit,r_value_re,r_value_neg_imag,L_CL,fitting_model)
            if output_full_analytical_impedance:
                Z_fit_full = v_analytical_H2_N2_impedance_Kulikovsky(sigma_e_fit,sigma_p_fit,L_CL,Cdl_cm3,frequencies_input)
                save_impedance_to_file(output_filename,frequencies_input,np.real(Z_fit_full),[-x for x in np.imag(Z_fit_full)])
            else:
                save_impedance_to_file(output_filename,frequencies_input_trunc,Z_re_fit,Z_neg_imag_fit)
    elif fitting_model=="H2_N2_graphical":
        se=Symbol('se')
        sp=Symbol('sp')
        
        # Fit the straight low-frequency part
        slope_LF, intercept_LF, r_value_lf, *tmp = stats.linregress(Z_re_input[:n_LF_points_to_use],Z_neg_imag_input[:n_LF_points_to_use])
        # Generate the straight line for LF
        x_LF = np.linspace(-intercept_LF/slope_LF,max(Z_re_input),2)
        y_LF = [intercept_LF+slope_LF*x for x in x_LF]
        
        LF_fit_slope = np.arctan(slope_LF)/np.pi*180
        if LF_fit_slope>87:
            print("Slope of the low-frequency fit is ",LF_fit_slope," degrees.")
        else:
            print("[Warning] Slope of the low-frequency fit is ",LF_fit_slope," degrees. The fitted results may be inaccurate.")
            
        slope_HF, intercept_HF, r_value_hf, *tmp =stats.linregress(Z_re_input_trunc,Z_neg_imag_input_trunc)
        # Generate the straight line for HF
        x_HF = np.linspace(-intercept_HF/slope_HF,max(Z_re_input_trunc),2)
        y_HF = [intercept_HF+slope_HF*x for x in x_HF]
        
        HF_fit_slope = np.arctan(slope_HF)/np.pi*180
        print("Slope of the high-frequency fit is ",HF_fit_slope," degrees.")
        
        if x_HF[0]<0:
            value = min(Z_re_input)
        else:
            value = x_HF[0]
        
        # Solve the system of equations for conductivities (see reference [2])
        sigma_pairs = solve((L_CL/(se+sp)-value/1000,L_CL/3/se+L_CL/3/sp-x_LF[0]/1000),se,sp)
        sigma_e_fit=float(sigma_pairs[0][0])
        sigma_p_fit=float(sigma_pairs[0][1])
        
        output_fitting_info(sigma_e_fit,sigma_p_fit,r_value_lf,r_value_hf,L_CL,fitting_model)
        Z_fit_full = v_analytical_H2_N2_impedance_Kulikovsky(sigma_e_fit,sigma_p_fit,L_CL,Cdl_cm3,frequencies_input)
        save_impedance_to_file(output_filename,frequencies_input,np.real(Z_fit_full),[-x for x in np.imag(Z_fit_full)])
    else:
        print("[Error] Unsupported fitting method.")
        sys.exit(1)
    
    # ------------ #
    # Nyquist plot #
    # ------------ #
    
    fig1 = plt.figure(1,facecolor='white', figsize=(500/dpi_val, 400/dpi_val), dpi=dpi_val)
    
    ax1 = plt.subplot2grid((1, 1), (0, 0))
    
    if fitting_model=="H2_O2" or fitting_model=="H2_N2":
        ax1.plot(Z_re_input_trunc,Z_neg_imag_input_trunc,marker="o",fillstyle="none",color="k",linestyle="none",markersize=10,label='Original data', zorder=5)
        ax1.plot(Z_re_fit,Z_neg_imag_fit,"r--",linewidth="3",label='Fitted impedance', zorder=10)
    elif fitting_model=="H2_N2_graphical":
        ax1.plot(Z_re_input,Z_neg_imag_input,marker="o",fillstyle="none",color="k",linestyle="none",markersize=10,label='Original data', zorder=5)
        ax1.plot(np.real(Z_fit_full),[-x for x in np.imag(Z_fit_full)],"r--",linewidth="3",label='Fitted impedance', zorder=10)
        # Plot the linear high-frequency and low-frequency fits
        ax1.plot(x_HF,y_HF,"b-",linewidth="3",label="High-frequency fit",zorder=11)
        ax1.plot(x_LF,y_LF,"b-.",linewidth="3",label="Low-frequency fit",zorder=10)
        # Place markers where the linear fits intercept with the real axis
        ax1.scatter(x_LF[0],0,marker="o",color="b",s=100,linewidth=3,zorder=100)
        ax1.scatter(x_HF[0],0,marker="o",color="b",s=100,linewidth=3,zorder=100)
        
        pad=(min(x_LF)-min(x_HF))/4.0
        ax1.set_xlim(min(x_HF)-pad,max(x_LF)+pad)
        ax1.set_ylim(min(x_HF)-pad,max(x_LF)+pad)
    
    ax1.set_xlabel(r'Re(Z), mOhm*cm2')
    ax1.set_ylabel(r'-Im(Z), mOhm*cm2')
    
    ax1.set_aspect('equal')
    
    ax1.legend(loc='best')
    
    fig1.tight_layout()
    #plt.show()
    
    # --------- #
    # Bode plot #
    # --------- #
    
    fig2 = plt.figure(2,facecolor='white', figsize=(500/dpi_val, 400/dpi_val), dpi=dpi_val)
    
    ax2 = plt.subplot2grid((1, 1), (0, 0))
    
    if fitting_model=="H2_O2" or fitting_model=="H2_N2":
        ax2.plot(frequencies_input_trunc,Z_neg_imag_input_trunc,marker="o",fillstyle="none",color="k",linestyle="none",markersize=10,label='Original data', zorder=5)
        ax2.plot(frequencies_input_trunc,Z_neg_imag_fit,"r--",linewidth="3",label='Fitted impedance', zorder=10)
    elif fitting_model=="H2_N2_graphical":
        ax2.plot(frequencies_input,Z_neg_imag_input,marker="o",fillstyle="none",color="k",linestyle="none",markersize=10,label='Original data', zorder=5)
        ax2.plot(frequencies_input,[-x for x in np.imag(Z_fit_full)],"r--",linewidth="3",label='Fitted impedance', zorder=10)
        
    ax2.set_xlabel(r'Frequency_Hz')
    ax2.set_ylabel(r'-Im(Z), mOhm*cm2')
    
    if Bode_log_scale:
        ax2.set_xscale('log')
        locmaj = ticker.LogLocator(base=10,numticks=15) 
        ax2.xaxis.set_major_locator(locmaj)
        locmin = ticker.LogLocator(base=10.0,subs=(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),numticks=12)
        ax2.xaxis.set_minor_locator(locmin)
        ax2.xaxis.set_minor_formatter(ticker.NullFormatter())
    
    ax2.legend(loc='best')
    
    fig2.tight_layout(pad=1.5)
    if fitting_model=="H2_N2_graphical":
        plt.show()
    
    if fitting_model=="H2_O2" or fitting_model=="H2_N2":
        # ------------- #
        # Residual plot #
        # ------------- #
        
        fig3 = plt.figure(3,facecolor='white', figsize=(500/dpi_val, 400/dpi_val), dpi=dpi_val)
        colorbar_pad=0.05
        formatstring=r'{:.2e}'
        n_levels=100
        fit_color="w"
        cbar_orientation='vertical'
        
        ax3 = plt.subplot2grid((1, 1), (0, 0))
        
        # Get conductivity data (X and Y) from the fit
        tmp = np.array([x[0] for x in residual[2][0]])
        X = np.linspace(min(tmp),max(tmp),100)
        Y = np.linspace(min(residual[2][1][0]),max(residual[2][1][0]),100)
        min_x = 0.5*sigma_e_fit
        max_x = 1.5*sigma_e_fit
        min_y = 0.5*sigma_p_fit
        max_y = 1.5*sigma_p_fit
        ind_x = np.argwhere((X>min_x-1e-14) & (X<max_x+1e-14))
        ind_y = np.argwhere((Y>min_y-1e-14) & (Y<max_y+1e-14))
        ind_x = ind_x.reshape(len(ind_x))
        ind_y = ind_y.reshape(len(ind_y))
        X = X[ind_x]
        Y = Y[ind_y]
        
        # Get residual data (Z) from the fit
        Z = np.zeros([len(X),len(Y)])
        for i in range(0,len(X)):
            for j in range(0,len(Y)):
                if fitting_model=="H2_O2":
                    Z[i,j] = residual_analytical_H2_O2_impedance_Kulikovsky([X[i],Y[j]],Z_re_input_trunc,Z_neg_imag_input_trunc,L_CL,Cdl_cm3,frequencies_input_trunc,dc_current_density,Tafel_slope,residual_method)
                elif fitting_model=="H2_N2":
                    Z[i,j] = residual_analytical_H2_N2_impedance_Kulikovsky([X[i],Y[j]],Z_re_input_trunc,Z_neg_imag_input_trunc,L_CL,Cdl_cm3,frequencies_input_trunc,residual_method)
        
        # Define a log scale for the colormap
        logmin = round(math.log10(min(np.min(Z),residual[1])))
        logmax = round(math.log10(np.max(Z)))
        norm = colors.LogNorm(vmin=10**logmin, vmax=10**logmax)
        # Plot the residual on the conductivity grid
        surf = ax3.contourf(Y, X, Z, levels=np.logspace(logmin,logmax,n_levels), locator=ticker.LogLocator(), norm=norm, extend="both")
        # Show and annotate the fitted point
        ax3.plot(sigma_p_fit,sigma_e_fit,marker="o",color=fit_color,markersize=10)
        ax3.annotate(formatstring.format(residual[1]),xy=(sigma_p_fit*1.05,sigma_e_fit*1.05),color=fit_color)
        
        ax3.set_xlabel(r'Protonic conductivity, S/cm')
        ax3.set_ylabel(r'Electronic conductivity, S/cm')  
        ticks = 10**np.linspace(logmin,logmax,-logmin+logmax+1)
        cbar = fig3.colorbar(surf, ax=ax3, shrink=1, aspect=20, orientation=cbar_orientation, pad=colorbar_pad, ticks=ticks)
        cbar.ax.set_ylabel('Residual')
        for c in surf.collections:
            c.set_edgecolor('face')
        
        fig3.tight_layout(pad=1)
        plt.show()

if __name__ == "__main__":
   main(sys.argv[1:])