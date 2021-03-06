#!/bin/python

###################################################################################################
### Name: Read zeTorch Data ### Developer: Tucker McClanahan
### Purpose: Reads zeTorch Data from a run where all of the run files are in the same directory
### Use: python read_zeTorch_data.py data_directory/ starting_wavelength ending_wavelength
### Dependents: data_directory/ --> this is the directory with all of the spectrometer data files
###                                 from a single run
###             starting_wavelength, ending_wavelength --> this is the starting and ending 
###                         wavelengths of a single pulse where these wavelengths define the pulse
###                         to be analyzed by this script.  
### Notes: 
###################################################################################################

import numpy as np
from scipy import signal as sig
from scipy import stats as stats
from scipy.special import wofz
import sys
import os
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from lmfit import Model, Parameters
from zeTorch.Run import Run 

def main():
    """ Reading and Processing the Run Data
    
    This definition is where the data is parsed and filtered by other definitions. It ends with the
    plotting of the data by another set of definitions. 
    """
    #path = '/Users/tuckermcclanahan/Google_Drive/PhD/Eriks_Data/zeTorch/Thesis Data/1 Low Position/01 Argon Pure/' 
    #directories = np.array(['20170803', 
    #                         '20170804', 
    #                         '20170807', 
    #                         '20170809', 
    #                         '20170906', 
    #                         '20171009']) 
    #filename = 'Low_Position_Argon_Pure.pdf'
    #pp = PdfPages(filename)
    # 
    #for i in range(np.size(directories)):
    #    directory_path = path+directories[i]+'/'
    #    print directory_path
    #    title = directories[i]
    #    all_data = parse_input(directory_path)
    #    for dat in all_data:
    #        dat.calibrate_data('./suspect_calibration_data/CalibrationFile.txt')
    #    corrected_data = fitting_bb_data(all_data)
    #    boltz_ar(corrected_data, pp, title)
    #
    #pp.close()
    # 
    #path = '/Users/tuckermcclanahan/Google_Drive/PhD/Eriks_Data/zeTorch/Thesis Data/1 Low Position/02 Argon Hyd Torch/'
    #directories = np.array(['20170803',  
    #                         '20170804', 
    #                         '20170807', 
    #                         '20170809', 
    #                         '20171009']) 
    #filename = 'Low_Position_Argon_Hyd_Torch.pdf'
    #pp = PdfPages(filename)
    # 
    #for i in range(np.size(directories)):
    #    directory_path = path+directories[i]+'/'
    #    print directory_path
    #    title = directories[i]
    #    all_data = parse_input(directory_path)
    #    for dat in all_data:
    #        dat.calibrate_data('./suspect_calibration_data/CalibrationFile.txt')
    #    corrected_data = fitting_bb_data(all_data)
    #    boltz_ar(corrected_data, pp, title)
    #    boltz_H(corrected_data, pp, title)
    #
    #pp.close() 
    #
    #path = '/Users/tuckermcclanahan/Google_Drive/PhD/Eriks_Data/zeTorch/Thesis Data/1 Low Position/03 Powder NY/'
    #directories = np.array(['20170803',  
    #                         '20170804', 
    #                         '20171009',
    #                         '20171101']) 
    #filename = 'Low_Position_Powder_NY.pdf'
    #pp = PdfPages(filename)
    # 
    #for i in range(np.size(directories)):
    #    directory_path = path+directories[i]+'/'
    #    print directory_path
    #    title = directories[i]
    #    all_data = parse_input(directory_path)
    #    for dat in all_data:
    #        dat.calibrate_data('./suspect_calibration_data/CalibrationFile.txt')
    #    corrected_data = fitting_bb_data(all_data)
    #    boltz_ar(corrected_data, pp, title)
    #    boltz_H(corrected_data, pp, title)
    #
    #pp.close() 
    #
    #path = '/Users/tuckermcclanahan/Google_Drive/PhD/Eriks_Data/zeTorch/Thesis Data/1 Low Position/04 Powder Trin/'
    #directories = np.array(['20170807',  
    #                         '20170809']) 
    #filename = 'Low_Position_Powder_Trin.pdf'
    #pp = PdfPages(filename)
    # 
    #for i in range(np.size(directories)):
    #    directory_path = path+directories[i]+'/'
    #    print directory_path
    #    title = directories[i]
    #    all_data = parse_input(directory_path)
    #    for dat in all_data:
    #        dat.calibrate_data('./suspect_calibration_data/CalibrationFile.txt')
    #    corrected_data = fitting_bb_data(all_data)
    #    boltz_ar(corrected_data, pp, title)
    #    boltz_H(corrected_data, pp, title)
    #
    #pp.close() 
    
    path = '/Users/tuckermcclanahan/Google_Drive/PhD/Eriks_Data/zeTorch/Thesis Data/1 Low Position/05 Powder Silica/'
    directories = np.array(['20170906'])  
    filename = 'Low_Position_Powder_Silica.pdf'
    pp = PdfPages(filename)
     
    for i in range(np.size(directories)):
        directory_path = path+directories[i]+'/'
        print directory_path
        title = directories[i]
        all_data = parse_input(directory_path)
        for dat in all_data:
            dat.calibrate_data('./suspect_calibration_data/CalibrationFile.txt')
        corrected_data = fitting_bb_data(all_data)
        boltz_ar(corrected_data, pp, title)
    
    pp.close() 
    
    path = '/Users/tuckermcclanahan/Google_Drive/PhD/Eriks_Data/zeTorch/Thesis Data/2 Elevated/01 Argon Pure/'
    directories = np.array(['20171114',
                            '20171201',
                            '20180221',
                            '20180222',
                            '20180320'])
    
      
    filename = 'Elevated_Argon_Pure.pdf'
    pp = PdfPages(filename)
     
    for i in range(np.size(directories)):
        directory_path = path+directories[i]+'/'
        print directory_path
        title = directories[i]
        all_data = parse_input(directory_path)
        for dat in all_data:
            dat.calibrate_data('./suspect_calibration_data/CalibrationFile.txt')
        corrected_data = fitting_bb_data(all_data)
        boltz_ar(corrected_data, pp, title)
    
    pp.close() 
    
    path = '/Users/tuckermcclanahan/Google_Drive/PhD/Eriks_Data/zeTorch/Thesis Data/2 Elevated/02 Argon Hyd Torch/'
    directories = np.array(['20171113', 
                            '20171114', 
                            '20171201', 
                            '20180221', 
                            '20180222', 
                            '20180320']) 
      
    filename = 'Elevated_Argon_Hyd_Torch.pdf'
    pp = PdfPages(filename)
     
    for i in range(np.size(directories)):
        directory_path = path+directories[i]+'/'
        print directory_path
        title = directories[i]
        all_data = parse_input(directory_path)
        for dat in all_data:
            dat.calibrate_data('./suspect_calibration_data/CalibrationFile.txt')
        corrected_data = fitting_bb_data(all_data)
        boltz_ar(corrected_data, pp, title)
        boltz_H(corrected_data, pp, title)
    
    pp.close() 
    
    path = '/Users/tuckermcclanahan/Google_Drive/PhD/Eriks_Data/zeTorch/Thesis Data/2 Elevated/03 Powder NY/'
    directories = np.array(['20171113', 
                            '20171114', 
                            '20171201', 
                            '20180222']) 
      
    filename = 'Elevated_Powder_NY.pdf'
    pp = PdfPages(filename)
     
    for i in range(np.size(directories)):
        directory_path = path+directories[i]+'/'
        print directory_path
        title = directories[i]
        all_data = parse_input(directory_path)
        for dat in all_data:
            dat.calibrate_data('./suspect_calibration_data/CalibrationFile.txt')
        corrected_data = fitting_bb_data(all_data)
        boltz_ar(corrected_data, pp, title)
        boltz_H(corrected_data, pp, title)
    
    pp.close() 
    
    path = '/Users/tuckermcclanahan/Google_Drive/PhD/Eriks_Data/zeTorch/Thesis Data/2 Elevated/04 Powder Trin/'
    directories = np.array(['20180320']) 
      
    filename = 'Elevated_Powder_Trin.pdf'
    pp = PdfPages(filename)
     
    for i in range(np.size(directories)):
        directory_path = path+directories[i]+'/'
        print directory_path
        title = directories[i]
        all_data = parse_input(directory_path)
        for dat in all_data:
            dat.calibrate_data('./suspect_calibration_data/CalibrationFile.txt')
        corrected_data = fitting_bb_data(all_data)
        boltz_ar(corrected_data, pp, title)
        boltz_H(corrected_data, pp, title)
    
    pp.close() 
    #directory_path = sys.argv[1]    

    #### Parse the files in the directory
    #all_data = parse_input(directory_path)
    #for dat in all_data:
    #    #dat.calibrate_data('./suspect_calibration_data/Calibrated_Source_Spectral_Output.txt', './suspect_calibration_data/DH-3PlusCalLight-DeuteriumHalogen_HRD10391_13-38-32-533.txt')
    #    dat.calibrate_data('./suspect_calibration_data/CalibrationFile.txt')
    ##pw = 655
    ##refined_data = refine_data(all_data, pw-10, pw+10)

    ##corrected_data = correct_the_data(all_data)
    #corrected_data = fitting_bb_data(all_data)
    ##data = temp_from_H(corrected_data)
    ##data = erics_module(corrected_data)
    #data = boltz_ar(corrected_data)
    #data = boltz_H(corrected_data)
#    gaussian_peak_fitting(all_data, peak_wavelength)
#    voigt_peak_fitting(all_data, peak_wavelength)
#    plot_filter(all_data)
#    plot_data(data, 'Plot')
    #plot_refined_data(refined_data, 'Plot')
    
    #print('Total number of data files: %d' % np.size(all_data))
     
#    refined_data = refine_data(all_data, starting_wavelength, ending_wavelength)
    
#    print 'Total number of Refined data files: %d' % np.size(refined_data)
    
#    plot_refined_data(refined_data, directory_path)
#    plot_data(refined_data, directory_path)

def boltz_H(data, pp, plot_title ):

    y_point = []
    x_point = []
    for dat in data:
        lam = np.asarray(dat.calibrated_lam)
        raw_counts = dat.calibrated_counts - np.asarray(dat.bb)
        #raw_counts = dat.calibrated_counts
        counts = raw_counts/np.trapz(raw_counts, x=lam)/0.473 # units of 1/nm 
         
        peaks_of_interest = np.array([656., 485.])
        max_peak_height = []
        max_peak_lam = []
        for poi in peaks_of_interest:
            peak_indexes = [i for i in range(np.size(lam)) if lam[i]>(poi-5.) and
            lam[i]<(poi+5.)]
            peak_lam = lam[peak_indexes[0]:peak_indexes[-1]] 
            peak_counts = counts[peak_indexes[0]:peak_indexes[-1]] 
            
            p = Parameters()
            p.add_many(('sigma'     ,        5.0,    True,    0.0,    None,    None),
                       ('gamma'     ,       1.0,    True,    0.0,      None,    None),
                       ('amp'     ,        0.1,    True,    0.0,    None,    None),
                       ('lam0'     ,        poi,    True,    0.0,    None,    None))
            func = Model(voigt)
            result = func.fit(peak_counts, lam=peak_lam, params=p)
            #print result.fit_report()

            pi = [i for i in range(np.size(peak_lam)) if peak_lam[i]>(poi-2.) and
            peak_lam[i]<(poi+2.)]
            #max_peak_height.append(np.trapz(result.eval(lam=lam[peak_indexes[0]:peak_indexes[-1]]), x= lam[peak_indexes[0]:peak_indexes[-1]]))
            #max_peak_height.append(np.trapz(result.best_fit, x = peak_lam))
            max_peak_height.append(np.trapz(peak_counts[pi[0]:pi[-1]], x = peak_lam[pi[0]:pi[-1]]))
            #max_peak_height.append(np.max(result.best_fit))
            max_peak_lam.append(result.params['lam0'].value)
            
            #plt.figure()
            #plt.xlabel('Wavelength (nm)')
            #plt.ylabel('Intensity (1/nm)')
            #plt.plot(peak_lam,peak_counts, 'b*')
            #plt.plot(peak_lam, result.best_fit, 'k-')
            #plt.savefig(pp, format='pdf')

        max_peak_lam = np.asarray(max_peak_lam)        
        max_peak_height = np.asarray(max_peak_height)
        g_k = np.array([18., 32.])
        A = np.array([ 0.441, 0.0841])
        e_k = np.array([97492., 102824.])
        e_i = np.array([82259., 82259.])
        e = e_k-e_i
        y = np.log(max_peak_height*max_peak_lam/(g_k*A))
        
        #indexs = np.array([1, 3, 7, 8])
        #y_point = np.concatenate((y_point,np.asarray([y[i] for i in indexs] ))) 
        #x_point = np.concatenate((x_point,np.asarray([e[i] for i in indexs] ))) 
        y_point = np.concatenate((y_point,y )) 
        x_point = np.concatenate((x_point,e)) 

    e = np.asarray(x_point)
    y = np.asarray(y_point)
    #slope, intercept, r_value, p_value, std_err = stats.linregress(np.asarray(e_k),np.asarray(y))
    slope, intercept, r_value, p_value, std_err = stats.linregress(e[e.argsort()],y[e.argsort()])
    
    boltzmann = 0.69503476
    temp = -1./(slope*boltzmann)
    #y_line = slope*np.asarray(e_k)+intercept
    y_line = slope*np.asarray(e)+intercept
    
    plt.figure()
    #plt.plot(e_k,y, '*')
    #plt.plot(e_k, y_line, '-')
    #plt.title(' Calc Temp = %1.4E K' % (temp))
    plt.title(plot_title+' H Calc Temp = %1.4E K' % (temp))
    plt.xlabel('Delta E (1/cm)')
    plt.ylabel('$\ln((I\lambda)/(g_k*A))$')
    plt.plot(e,y, '*')
    plt.plot(e, y_line, '-')
    plt.savefig(pp, format='pdf')


def boltz_ar(data, pp, plot_title):

    y_point = []
    x_point = []
    for dat in data:
        lam = np.asarray(dat.calibrated_lam)
        raw_counts = dat.calibrated_counts - np.asarray(dat.bb)
        #raw_counts = dat.calibrated_counts
        counts = raw_counts/np.trapz(raw_counts, x=lam)/0.473 # units of 1/nm 
         
        peaks_of_interest = np.array([810.035, 762.182, 749.397, 737.032, 705.584, 695.367, 839.292, 601.445, 910.828]) 
        max_peak_height = []
        max_peak_lam = []
        for poi in peaks_of_interest:
            peak_indexes = [i for i in range(np.size(lam)) if lam[i]>(poi-5.) and
            lam[i]<(poi+5.)]
            peak_lam = lam[peak_indexes[0]:peak_indexes[-1]] 
            peak_counts = counts[peak_indexes[0]:peak_indexes[-1]] 
            
            p = Parameters()
            p.add_many(('sigma'     ,        5.0,    True,    0.0,    None,    None),
                       ('gamma'     ,       1.0,    True,    0.0,      None,    None),
                       ('amp'     ,        0.1,    True,    0.0,    None,    None),
                       ('lam0'     ,        poi,    True,    0.0,    None,    None))
            func = Model(voigt)
            result = func.fit(peak_counts, lam=peak_lam, params=p)
            #print result.fit_report()

            pi = [i for i in range(np.size(peak_lam)) if peak_lam[i]>(poi-2.) and
            peak_lam[i]<(poi+2.)]
            #max_peak_height.append(np.trapz(result.eval(lam=lam[peak_indexes[0]:peak_indexes[-1]]), x= lam[peak_indexes[0]:peak_indexes[-1]]))
            #max_peak_height.append(np.trapz(result.best_fit, x = peak_lam))
            max_peak_height.append(np.trapz(peak_counts[pi[0]:pi[-1]], x = peak_lam[pi[0]:pi[-1]]))
            #max_peak_height.append(np.max(result.best_fit))
            max_peak_lam.append(result.params['lam0'].value)
            
            #plt.figure()
            #plt.xlabel('Wavelength (nm)')
            #plt.ylabel('Intensity (1/nm)')
            #plt.plot(peak_lam,peak_counts, 'b*')
            #plt.plot(peak_lam, result.best_fit, 'k-')
            #plt.savefig(pp, format='pdf')

        max_peak_lam = np.asarray(max_peak_lam)        
        max_peak_height = np.asarray(max_peak_height)
        g_k = np.array([7.0,     5.,      1.0,     5.,      5.,       3.,     5.,      9.,      3.]) # no units
        A = np.array([ 0.366,    0.274,   0.472,   0.087,   0.0395,  0.067,  0.244,    0.0246,   0.212]) # 1E8 1/s
        e_k = np.array([105463., 106238., 108723., 107290., 107290., 107496., 107290., 122036., 104102.]) # 1/cm
        e_i = np.array([93144.,  93144.,  95400.,  93751.,  93144.,  93144.,  95400.,  105463., 93144.])
        y = np.log(max_peak_height*max_peak_lam/(g_k*A))
        e = e_k-e_i
        
        indexs = np.array([1, 3, 7, 8])
        y_point = np.concatenate((y_point,np.asarray([y[i] for i in indexs] ))) 
        x_point = np.concatenate((x_point,np.asarray([e[i] for i in indexs] ))) 

    e = np.asarray(x_point)
    y = np.asarray(y_point)
    #slope, intercept, r_value, p_value, std_err = stats.linregress(np.asarray(e_k),np.asarray(y))
    slope, intercept, r_value, p_value, std_err = stats.linregress(e[e.argsort()],y[e.argsort()])
    
    boltzmann = 0.69503476
    temp = -1./(slope*boltzmann)
    #y_line = slope*np.asarray(e_k)+intercept
    y_line = slope*np.asarray(e)+intercept
    
    plt.figure()
    #plt.plot(e_k,y, '*')
    #plt.plot(e_k, y_line, '-')
    plt.title(plot_title+' Ar Calc Temp = %1.4E K' % (temp))
    plt.xlabel('Delta E (1/cm)')
    plt.ylabel('$\ln((I\lambda)/(g_k*A))$')
    plt.plot(e,y, '*')
    plt.plot(e, y_line, '-')
    plt.savefig(pp, format='pdf')

def voigt(lam, sigma, gamma, amp, lam0):
    #sigma = alpha/np.sqrt(2*np.log(2))
    z = ((lam-lam0)+gamma*1j)/(sigma*np.sqrt(2))
    return np.real(wofz(z))/(sigma*np.sqrt(2.*np.pi))*amp

def gaussian_peak_fitting(all_data, peak_wavelength):
    all_data = fitting_bb_data(all_data)
    for dat in all_data:
        counts = np.asarray(dat.corrected_counts)
        lam = np.asarray(dat.wavelengths)
        peak_indexes = [i for i in range(np.size(lam)) if lam[i]>(peak_wavelength-10.) and
        lam[i]<(peak_wavelength+10.)]
        
### Before BB substraction
        peak_lam = lam[peak_indexes[0]:peak_indexes[-1]] 
        peak_counts = counts[peak_indexes[0]:peak_indexes[-1]] 
        p = Parameters()
         #          (Name   ,        Value,    Vary,    Min,     Max,    Expr)
        p.add_many(('amp'     ,        15000,    True,    None,    None,    None),
                   ('cen'     ,        peak_wavelength,    True,    None,    None,    None),
                   ('wid'     ,        3.0,    True,    None,    None,    None),
                   ('scale',           0.0, True, None, None, None))

        func = Model(gaussian)
        result = func.fit(peak_counts, x=peak_lam, params=p)
        print(result.fit_report())

### After BB subtraction
        peak_bb_sub_counts = peak_counts - np.asarray(dat.bb[peak_indexes[0]:peak_indexes[-1]])
        p = Parameters()
         #          (Name   ,        Value,    Vary,    Min,     Max,    Expr)
        p.add_many(('amp'     ,        14000.,    True,    None,    None,    None),
                   ('cen'     ,        peak_wavelength,    True,    None,    None,    None),
                   ('wid'     ,        3.0,    True,    None,    None,    None),
                   ('scale',           0.0, True, None, None, None))

        func = Model(gaussian)
        result_bb = func.fit(peak_bb_sub_counts, x=peak_lam, params=p)
        print("With Black Body Subtraction")
        print(result_bb.fit_report())

        
 
        pp = PdfPages('Gaussian_Fit.pdf')
        plt.figure()
        #plt.plot(lam,counts, 'k-', label='Raw Data')
        plt.plot(peak_lam,peak_counts, 'ro', label='Peak Data')
        plt.plot(peak_lam, result.best_fit, 'b-', label='Best Fit Peak Data')
        plt.plot(peak_lam,peak_bb_sub_counts, 'ko', label='Peak BB Sub Data')
        plt.plot(peak_lam, result_bb.best_fit, 'm-', label='Best Fit Peak BB Sub Data')
        plt.legend()
        plt.savefig(pp, format='pdf')
        pp.close()
    
def gaussian(x, amp, cen, wid, scale):
    return amp*np.exp(-(x-cen)**2/wid)+scale


def fitting_bb_data(all_data):
    """ Fits the black body portion of the spectrum for each datafile.

    This definition goes through all of the Run classes for a given dataset and filters out the
    peaks to fit a black body curve to the filtered data. The LMFIT fitting routine is used as a wrapper for the SCIPY optmize tools to fit the Planck's Black Body curve. This function feeds in initial guesses for the parameters and returns a best fitted parameters for the curve.  Keyword Arguments:
    all_data -- List of Run classes

    """

    #cal_data = np.loadtxt(os.getcwd()+'/suspect_calibration_data/CalibrationFile.txt')
    filtered_data = []
    #count = 0
    for dat in all_data:
        counts = np.asarray(dat.calibrated_counts)
        lam = np.asarray(dat.calibrated_lam)
        bb = sig.medfilt(counts,kernel_size=81)

        p = Parameters()
         #          (Name   ,        Value,    Vary,    Min,     Max,    Expr)
        p.add_many(('T'     ,        5000.,    True,    None,    None,    None),
                   ('scale'     ,        1E-11,    True,    None,    None,    None),
                   ('shift'     ,        0.0,    False,    None,    None,    None))

        func = Model(pbb)
        result = func.fit(bb, lam=lam, params=p)
        #print(dat.filename)
        #print(result.fit_report())
        
        dat.bb = bb
        dat.temp = result.params['T'].value
        dat.temp_err = result.params['T'].stderr
        dat.aic = result.aic
        filtered_data.append(dat) 

       # pp = PdfPages('Filter_Test.pdf')
       # plt.figure()
       # plt.plot(lam,counts, 'k-', label='Raw Data')
       # plt.plot(lam[1::], corrected_counts, 'm-', label='Corrected Data')
       # ##plt.plot(cal_lam[1::], deconc_counts[1], 'm-', label='Deconc Data')
       # ##plt.plot(cal_lam[1::], rebinned_counts, 'c--', label='Deconc Data')
       # ##plt.plot(lam[1::], bb, 'r-', label='Filtered Data')
       # plt.plot(lam[1::], result.best_fit, 'b-', label='Fit Filtered Data')
       # plt.legend()
       # plt.savefig(pp, format='pdf')
      
       # plt.figure()
       ## plt.plot(lamp_det_lam[1::], cal_calculated, 'ko', label='Calculated') 
       # plt.plot(cal_lam, cal_counts, 'b*', label='From File') 
       # plt.legend()
       # plt.savefig(pp, format='pdf')
    #pp.close() 

    return filtered_data

def rebin(x1, y1, x2):
    # Newly binned data
    y2 = []
    for i in range(x2.size-1):
        bind = 0.0
        count = 0.0
        for t in range(x1.size):
            if x1[t]>x2[i] and x1[t]<x2[i+1]:
                count += 1.0
                bind += y1[t]
        if count == 0 and i == 0:
            y2.append(0.0)
        elif count == 0 and i >0 :
            y2.append(y2[-1])
        else:
            y2.append(bind/count)
    return np.asarray(y2)

def pbb(lam, T, scale, shift):
    h = 6.626070040E-34 #Js
    c = 299792458. #m/s
    k = 1.38064852E-23 #J/K
    lamm = lam*1E-9 #m
    return ((2.*h*c**2/(lamm+shift)**5)*(1./(np.exp(h*c/((lamm+shift)*k*T))-1.)))*scale



def plot_filter(all_data):
    pp = PdfPages('Filter_Test.pdf')
    savgol_filtered = sig.savgol_filter(all_data[0].counts, 51, 2, mode='nearest') 
    med_filtered = sig.medfilt(all_data[0].counts, kernel_size=81) 
    #f = open('out_data.txt', 'w')
    #f.write('### Test Erik data with BB subtracted\n')
    #f.write('### Wavelength   Counts_Original   Counts_BB_corrected\n')
    #f.write('###   (nm)       \n')
    #for i in range(np.size(all_data[0].wavelengths)):
    #    f.write('%1.4E  %1.4E  %1.4E \n' % (all_data[0].wavelengths[i], all_data[0].counts[i],
    #    all_data[0].counts[i]-med_filtered[i]))
    plt.figure()
    plt.plot(all_data[0].wavelengths, all_data[0].counts, label='Original')
    plt.plot(all_data[0].wavelengths, savgol_filtered, label='Savgol Filtered')
    plt.plot(all_data[0].wavelengths, med_filtered, label='Med Filtered')
    plt.legend()
    plt.savefig(pp, format='pdf')
    pp.close()

def plot_data(all_data, path):
    """ Making a heatmap plot of all the data in a run along with individual plots of each timestep.

    Keyword arguments:
    all_data -- List of Run classes
    path -- String containing the path to the data directory
    
    File Outputs:
    path+'_All_Data.pdf' -- PDF of the heatmap plot and the individual plots of each timestep.
    """
     
    path = path.replace('/','_')
    pp = PdfPages(path+'_All_Data.pdf')
    times = []
    lam = np.asarray(all_data[0].corrected_lam)
     
    L, T = np.meshgrid(lam, times)

    data = np.zeros([np.size(times), np.size(lam)])
    temps = []
    temps_err = []
    h_boltz_temps = []

    count = 0
    for dat in all_data:
        if count == 0:
            first_time = dat.time
            times.append(0.0)
            count += 1
        else:
            times.append(dat.time-first_time)
        temps.append(dat.temp)
        temps_err.append(dat.temp_err)
        h_boltz_temps.append(dat.H_boltz_temp)
    #    for l in range(np.size(lam)):
    #        data[t,l] = abs(all_data[t].corrected_counts[l])
    times = np.asarray(times)/1000.
    #plt.figure()
    #plt.title('Plot of All Data')
    #plt.xlabel('Wavelength (nm)')
    #plt.ylabel('Time (AU)')
    #plt.pcolormesh(T, L, data, rasterized=True, cmap='hot')
    #plt.colorbar(label='Counts')
    #plt.savefig(pp, format='pdf')

    #plt.figure()
    #plt.title('Plot of All Data')
    #plt.xlabel('Wavelength (nm)')
    #plt.ylabel('Time (AU)')
    #plt.zlabel('Counts')
    #Axes3D.plot_surface(T,L, data)
    #plt.savefig(pp, format='pdf')
    #plt.show()

### Make the time plot of temps
    plt.figure()
    plt.title('Time vs. Temperature from BB')
    plt.xlabel('Time (sec)')
    plt.ylabel('Temperature (K)')
    plt.errorbar(times,temps, yerr=temps_err, fmt='bo')
    plt.savefig(pp, format='pdf')

    plt.figure()
    plt.title('Time vs. Temperature from H Peaks')
    plt.xlabel('Time (sec)')
    plt.ylabel('Temperature (K)')
    plt.ylim((0.0, 10000))
    plt.plot(times,h_boltz_temps, 'bo')
    plt.savefig(pp, format='pdf')

    pp.close()


def plot_refined_data(refined_data, path):
    """ Plots the filtered data."""
    path = path.replace('/','_')
    pp = PdfPages(path+'.pdf')
    
    peak_heights = []
    peak_areas = []
    peak_fwhms = []

    for run in refined_data:
        peak_heights.append(run.max)
        peak_areas.append(run.area)
    #    peak_fwhms.append(run.fwhm)
        
### Make the time plot of peak_heights
    plt.figure()
    plt.title('Time vs. Peak Height')
    plt.xlabel('Time (au)')
    plt.ylabel('Peak Height (counts)')
    plt.plot(range(np.size(peak_heights)),peak_heights, 'ko')
    plt.savefig(pp, format='pdf')
 
### Make the time plot of peak_areas
    plt.figure()
    plt.title('Time vs. Peak Areas')
    plt.xlabel('Time (au)')
    plt.ylabel('Peak Areas (nm*counts)')
    plt.plot(range(np.size(peak_areas)),peak_areas, 'bo')
    plt.savefig(pp, format='pdf')


### Make the time plot of peak_fwhms
    #plt.figure()
    #plt.title('Time vs. Peak FWHM')
    #plt.xlabel('Time (au)')
    #plt.ylabel('Peak FWHM (counts)')
    #plt.plot(range(np.size(peak_fwhms)),peak_fwhms, 'bo')
    #plt.savefig(pp, format='pdf')
   
    pp.close() 

def refine_data(all_data, p_s, p_e):
    """ Filters down the data based on some selection criteria. 

    Keyword arguments:
    all_data -- List of Run classes
    p_s -- integer of the starting wavelength of a peak
    p_e -- integer of the ending wavelength of a peak

    Returns:
    list of Run classes of the runs that meet the selection criteria
    """
    data = []
    indexes = [i for i in range(np.size(all_data[0].wavelengths)) if all_data[0].wavelengths[i]>p_s and all_data[0].wavelengths[i]<p_e]
    for run in all_data:
        ## Find the peak around 655 +- 10
        peak_lam = run.wavelengths[indexes[0]:indexes[-1]]
        peak_data = run.counts[indexes[0]:indexes[-1]]
        peak_height = np.amax(peak_data)
        if peak_height < 500:
            continue
        
        run.max = peak_height
        run.area = np.trapz(peak_data,peak_lam)
        #run.fwhm = FWHM(peak_lam, peak_data)
	
	data.append(run)

    return data


def FWHM(x,y):
    """ Calculates the FWHM of a given pulse

    Keyword Arguments:
    x -- List of floats 
    y -- List of floats
    """
    half_max = max(y)/2.
    d = np.sign(half_max - np.array(y[0:-1])) - np.sign(half_max - np.array(y[1:]))
    left_idx = np.argwhere(d > 0)[0,0]
    right_idx = np.argwhere(d<0)[0,-1]
    return x[right_idx]-x[left_idx] 


        


    
     
    
def parse_input(dirpath):
    """ Parse the datafiles in the data directory.
    
    This definition goes through each datafile in the data directory. It goes line by line in each
    datafile looking for specific keywords and splits the line and assigns certain values in the
    line to an attribute in the Run() class. When every file in the data directory has been gone
    through, the definintion returns a list of Run classes where each Run class corresponds to a
    specific data file within the data directory. 
     
    Keyword arguments:
    dirpath -- string of the path to data directory
    
    Returns:
    List of Run classes
    """

    if not dirpath.endswith('/'):
        dirpath = dirpath+'/'
    filenames = [f for f in sorted(os.listdir(dirpath)) if os.path.isfile(os.path.join(dirpath, f))]
    
    data = []

    for filename in filenames:
        if filename.startswith('.'):
            continue
        run = Run()
        run.load_file(os.path.join(dirpath,filename))
        data.append(run)

    return data

def run_checks():
    """ Making sure the user inputs the correct inputs


    Making sure that the user inputs the correct number of expected inputs on the commmand line and
    that those inputs are of the correct type and can be read by PYTHON. 
    """
    if np.size(sys.argv)<4:
        raise Exception("Please give python read_erics_data.py data_directory_name starting_wavelength ending_wavelength)")
    dir_path = sys.argv[1]
    try:
        sw = float(sys.argv[2])
    except ValueError:
        print "Please input a number as the third entry"
    try:
        ew = float(sys.argv[3])
    except ValueError:
        print "Please input a number as the fourth entry"
    if not os.path.isdir(dir_path): 
        raise Exception("Directory not found!")
        if not os.access(dir_path, os.R_OK):
            raise Exception("Directory is unreadable")

    return dir_path, sw, ew
    



main()
