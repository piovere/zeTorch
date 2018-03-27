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

def main():
    """ Reading and Processing the Run Data
    
    This definition is where the data is parsed and filtered by other definitions. It ends with the
    plotting of the data by another set of definitions. 
    """
    path = '/Users/tuckermcclanahan/Google_Drive/PhD/Eriks_Data/zeTorch/Thesis Data/1 Low Position/01 Argon Pure/' 
    directories = np.array(['20170803']) 
    #                         '20170804', 
    #                         '20170807', 
    #                         '20170809', 
    #                         '20170906', 
    #                         '20171009']) 
    filename = 'Low_Position_Argon_Pure.pdf'
    pp = PdfPages(filename)
     
    for i in range(np.size(directories)):
        directory_path = path+directories[i]+'/'
        print directory_path
        title = directories[i]
        all_data = parse_input(directory_path)
        corrected_data = fitting_bb_data(all_data)
        boltz_ar(corrected_data, pp, title)
    
    pp.close()
     
    path = '/Users/tuckermcclanahan/Google_Drive/PhD/Eriks_Data/zeTorch/Thesis Data/1 Low Position/02 Argon Hyd Torch/'
    directories = np.array(['20170803',  
                             '20170804', 
                             '20170807', 
                             '20170809', 
                             '20171009']) 
    filename = 'Low_Position_Argon_Hyd_Torch.pdf'
    pp = PdfPages(filename)
     
    for i in range(np.size(directories)):
        directory_path = path+directories[i]+'/'
        print directory_path
        title = directories[i]
        all_data = parse_input(directory_path)
        corrected_data = fitting_bb_data(all_data)
        boltz_ar(corrected_data, pp, title)
    
    pp.close() 
    
    path = '/Users/tuckermcclanahan/Google_Drive/PhD/Eriks_Data/zeTorch/Thesis Data/1 Low Position/03 Powder NY/'
    directories = np.array(['20170803',  
                             '20170804', 
                             '20171009',
                             '20171101']) 
    filename = 'Low_Position_Powder_NY.pdf'
    pp = PdfPages(filename)
     
    for i in range(np.size(directories)):
        directory_path = path+directories[i]+'/'
        print directory_path
        title = directories[i]
        all_data = parse_input(directory_path)
        corrected_data = fitting_bb_data(all_data)
        boltz_ar(corrected_data, pp, title)
    
    pp.close() 
    
    path = '/Users/tuckermcclanahan/Google_Drive/PhD/Eriks_Data/zeTorch/Thesis Data/1 Low Position/04 Powder Trin/'
    directories = np.array(['20170807',  
                             '20170809']) 
    filename = 'Low_Position_Powder_Trin.pdf'
    pp = PdfPages(filename)
     
    for i in range(np.size(directories)):
        directory_path = path+directories[i]+'/'
        print directory_path
        title = directories[i]
        all_data = parse_input(directory_path)
        corrected_data = fitting_bb_data(all_data)
        boltz_ar(corrected_data, pp, title)
    
    pp.close() 
    #directory_path = sys.argv[1]    

    ### Parse the files in the directory
    #all_data = parse_input(directory_path)

    #pw = 655
    #refined_data = refine_data(all_data, pw-10, pw+10)

    #corrected_data = correct_the_data(all_data)
    #corrected_data = fitting_bb_data(all_data)
    #data = temp_from_H(corrected_data)
    #data = erics_module(corrected_data)
    #data = boltz_ar(corrected_data)
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

def boltz_ar(data, pp, plot_title):

    y_point = []
    x_point = []
    for dat in data:
        lam = np.asarray(dat.corrected_lam)
        raw_counts = np.asarray(dat.corrected_counts)- np.asarray(dat.bb)
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
    plt.title(plot_title+' Calc Temp = %1.4E K' % (temp))
    plt.xlabel('Delta E (1/cm)')
    plt.ylabel('$\ln((I\lambda)/(g_k*A))$')
    plt.plot(e,y, '*')
    plt.plot(e, y_line, '-')
    plt.savefig(pp, format='pdf')
    

def erics_module(data):

    y_point = []
    x_point = []
    for dat in data:
        lam = np.asarray(dat.corrected_lam)
        counts = np.asarray(dat.corrected_counts)- np.asarray(dat.bb)
        
        pp = PdfPages('Argon_Temp_Calc.pdf')
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
                       ('amp'     ,        3E4,    True,    0.0,    None,    None),
                       ('lam0'     ,        poi,    True,    0.0,    None,    None))
            func = Model(voigt)
            result = func.fit(peak_counts, lam=peak_lam, params=p)
            print result.fit_report()

            peak_indexes = [i for i in range(np.size(lam)) if lam[i]>(poi-2.) and
            lam[i]<(poi+2.)]
            #max_peak_height.append(np.trapz(result.eval(lam=lam[peak_indexes[0]:peak_indexes[-1]]), x= lam[peak_indexes[0]:peak_indexes[-1]]))
            #max_peak_height.append(np.trapz(result.best_fit, x = peak_lam))
            max_peak_height.append(np.trapz(peak_counts, x = peak_lam))
            max_peak_lam.append(result.params['lam0'].value)
            plt.figure()
            plt.xlabel('Wavelength (nm)')
            plt.ylabel('Counts')
            plt.plot(peak_lam,peak_counts, 'b*')
            plt.plot(peak_lam, result.best_fit)
            plt.savefig(pp, format='pdf')

        max_peak_lam = np.asarray(max_peak_lam)        
        max_peak_height = np.asarray(max_peak_height)
        
        ### Source: Atomic Transition Probabilities Na- Ca NSRDS-NBS 22 Volume II
        #g_k = np.array([7.0, 5., 1.0, 5., 5., 3., 5., 9., 3.]) # no units
        #A = np.array([ 0.366,    0.274,   0.472,   0.087,   0.0395,  0.067,  0.244,    0.0246,   0.212]) # 1E8 1/s
        #e_k = np.array([105463., 106238., 108723., 107290., 107290., 107496., 107290., 122036., 104102.]) # 1/cm
        
        ###  
        g_k = np.array([7.0, 5., 1.0, 5., 5., 3., 5., 9., 3.]) # no units
        A = np.array([ 0.366,    0.274,   0.472,   0.087,   0.0395,  0.067,  0.244,    0.0246,   0.212]) # 1E8 1/s
        e_k = np.array([105463., 106238., 108723., 107290., 107290., 107496., 107290., 122036., 104102.]) # 1/cm

        y = np.log(max_peak_height*max_peak_lam/(g_k*A))
        

        for i in range(np.size(peaks_of_interest)):
            count = 0
            for j in range(np.size(peaks_of_interest)):
                if i==j or j<i or i==6 or j==6 or i==7 or j==7:
                    continue
                print i, j
                y_point.append(y[i]-y[j])
                x_point.append(-1.*(e_k[i]-e_k[j]))
        
    slope, intercept, r_value, p_value, std_err = stats.linregress(np.asarray(x_point),np.asarray(y_point))
    
    boltzmann = 0.69503476
    temp = -1./(slope*boltzmann)
    print temp
     
    y = slope*np.asarray(x_point)+intercept
    
    plt.figure()
    plt.ylim([-1,1])
    plt.plot(x_point,y_point, '*')
    plt.plot(x_point, y, '-')
    plt.savefig(pp, format='pdf')
    
    pp.close()

def temp_from_H(data):
    
    boltz_temps = [] 
    for dat in data:
        lam = np.asarray(dat.corrected_lam)
        counts = np.asarray(dat.corrected_counts)
       
        pp = PdfPages('Voigt_Fit.pdf')
        peaks_of_interest = np.array([656., 486.])
        max_peak_height = []
        max_peak_lam = []
        for poi in peaks_of_interest:
            peak_indexes = [i for i in range(np.size(lam)) if lam[i]>(poi-10.) and
            lam[i]<(poi+10.)]
            peak_lam = lam[peak_indexes[0]:peak_indexes[-1]] 
            peak_counts = counts[peak_indexes[0]:peak_indexes[-1]] 
            
            p = Parameters()
            p.add_many(('sigma'     ,        5.0,    True,    0.0,    None,    None),
                       ('gamma'     ,       1.0,    True,    0.0,      None,    None),
                       ('amp'     ,        3E4,    True,    0.0,    None,    None),
                       ('lam0'     ,        poi,    True,    0.0,    None,    None))
            func = Model(voigt)
            result = func.fit(peak_counts, lam=peak_lam, params=p)
            print result.fit_report()

            max_peak_height.append(np.max(result.best_fit))
            max_peak_lam.append(peak_lam[np.argmax(result.best_fit)])
            plt.figure()
            plt.xlabel('Wavelength (nm)')
            plt.ylabel('Counts')
            plt.plot(peak_lam,peak_counts, 'b*')
            plt.plot(peak_lam, result.best_fit)
            plt.savefig(pp, format='pdf')

        pp.close()
        max_peak_lam = np.asarray(max_peak_lam)        
        max_peak_height = np.asarray(max_peak_height)
        g_weights = np.array([18., 32.])
        trans_prob = np.array([0.441, 0.08419])
        energy_trans = np.array([-2.4208E-19, -1.3617E-19])

        boltz_const = 1.38E-23

        fd = energy_trans*max_peak_height/(g_weights*trans_prob)
        dat.H_boltz_temp = -(energy_trans[0]-energy_trans[1])/np.log(fd[0]/fd[1])/boltz_const
        
    return data 
             
        

def voigt(lam, sigma, gamma, amp, lam0):
    #sigma = alpha/np.sqrt(2*np.log(2))
    z = ((lam-lam0)+gamma*1j)/(sigma*np.sqrt(2))
    return np.real(wofz(z))/(sigma*np.sqrt(2.*np.pi))*amp

def correct_the_data(all_data):
    cal_data = np.loadtxt(os.getcwd()+'/suspect_calibration_data/CalibrationFile.txt')
    cal_lam = cal_data[:,0].astype(float)
    cal_counts = cal_data[:,1].astype(float)
    dat = all_data[0]
    lam = np.asarray(dat.wavelengths)
    rebinned_cal_counts = rebin(cal_lam, cal_counts, lam)
    corrected_data = []
    for dat in all_data:
        counts = np.asarray(dat.counts)
        dat.corrected_counts = counts[1::]*rebinned_cal_counts
        dat.corrected_lam = lam[1::]
        corrected_data.append(dat) 
    return corrected_data 

def voigt_peak_fitting(all_data, peak_wavelength):
    all_data = fitting_bb_data(all_data)
    for dat in all_data:
        counts = np.asarray(dat.corrected_counts)
        lam = np.asarray(dat.wavelengths)
        peak_indexes = [i for i in range(np.size(lam)) if lam[i]>(peak_wavelength-10.) and
        lam[i]<(peak_wavelength+10.)]
        peak_lam = lam[peak_indexes[0]:peak_indexes[-1]] 
        peak_counts = counts[peak_indexes[0]:peak_indexes[-1]] 
        
### After BB subtraction
        peak_bb_sub_counts = peak_counts - np.asarray(dat.bb[peak_indexes[0]:peak_indexes[-1]])
        p = Parameters()
         #          (Name   ,        Value,    Vary,    Min,     Max,    Expr)
        p.add_many(('alpha'     ,        3.0,    True,    0.0,    None,    None),
                   ('gamma'     ,       3.0,    True,    0.0,      None,    None),
                   ('amp'     ,        70000.,    True,    0.0,    None,    None),
                   ('lam0'     ,        peak_wavelength,    True,    None,    None,    None))

        func = Model(voigt)
        result_bb = func.fit(peak_counts, lam= peak_lam, params=p)
        print( "With Black Body Subtraction")
        print(result_bb.fit_report())

        
 
        pp = PdfPages('Voigt_Fit.pdf')
        plt.figure()
        plt.plot(peak_lam,peak_bb_sub_counts, 'ko', label='Peak BB Sub Data')
        plt.plot(peak_lam, result_bb.best_fit, 'm-', label='Best Fit Peak BB Sub Data')
        plt.legend()
        plt.savefig(pp, format='pdf')
        pp.close()


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
    cal_data=np.loadtxt(os.getcwd()+'/suspect_calibration_data/Calibrated_Source_Spectral_Output.txt',delimiter=',',skiprows=1)
    lamp_det = Run()
    lamp_det.load_file(os.getcwd()+'/suspect_calibration_data/DH-3PlusCalLight-DeuteriumHalogen_HRD10391_13-38-32-533.txt')
    lamp_det_lam = np.asarray(lamp_det.wavelengths)
    lamp_det_counts = np.asarray(lamp_det.counts)
    cal_lam = cal_data[:,0].astype(float)
    cal_counts = cal_data[:,1].astype(float)
    #source_lam = source_data[:,0].astype(float)
    #source_counts = source_data[:,1].astype(float)
    
    filtered_data = []
    #count = 0
    for dat in all_data:
        #count += 1
        #print(count)
        counts = np.asarray(dat.counts)
        lam = np.asarray(dat.wavelengths)
        rebinned_cal_counts = rebin(cal_lam, cal_counts, lam)
        #rebinned_counts = rebin(lam, counts, cal_lam)
        rebinned_lamp_counts = rebin(lamp_det_lam, lamp_det_counts, lam)
        #fft_lamp_det_counts = np.fft.fft(lamp_det_counts[1::])
        #fft_lamp_source_counts = np.fft.fft(rebinned_source_counts)
        #cal_calculated = np.fft.ifft(fft_lamp_det_counts/fft_lamp_source_counts)
        corrected_counts = counts[1::]*rebinned_cal_counts/rebinned_lamp_counts
        #deconc_counts = sig.deconvolve(rebinned_counts, cal_counts[1::])
        bb = sig.medfilt(corrected_counts,kernel_size=81)

        p = Parameters()
         #          (Name   ,        Value,    Vary,    Min,     Max,    Expr)
        p.add_many(('T'     ,        5000.,    True,    None,    None,    None),
                   ('scale'     ,        1E-11,    True,    None,    None,    None),
                   ('shift'     ,        0.0,    False,    None,    None,    None))

        func = Model(pbb)
        result = func.fit(bb, lam=lam[1::], params=p)
        #print(dat.filename)
        #print(result.fit_report())
        
        dat.bb = bb
        dat.corrected_lam = lam[1::]
        dat.corrected_counts = corrected_counts
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


class Run(object):
    """ Class for one run datafile's contents.
    
    This class contains every piece of information for a single data file.

    Attributes:
    self.filename -- String of the filename
    self.time -- Float of the time data in milliseconds 
    self.name -- String of name within the file
    self.user -- String of the username
    self.spectrometer -- String of the spectrometer name
    self.trigger_mode -- Integer of the triggermode of the spectrometer
    self.integration_time -- Float of the integration time
    self.scans_to_average -- Integer of the scans to average value
    self.electric_dark_correction_enabled -- String
    self.nonlinearity_correction_enabled -- String
    self.boxcar_width -- Integer of the width of the boxcar
    self.xaxis_mode -- String
    self.number_of_pixels -- Integer
    self.wavelengths -- List of floats of the wavelengths in units of nm
    self.counts -- List of floats of the counts 
    
    self.add_data(self,line) -- Definition to add wavelength and count data to their respective
    lists
    self.load_file(self, filename) -- Loads the contents of the file and assigns it to the
    respective attribute 
    """   


    
    def __init__(self):
        self.filename = ''
        self.time = 0.0
        self.name = ''
        self.user = ''
        self.spectrometer = ''
        self.trigger_mode = 0
        self.integration_time = 0.0
        self.scans_to_average = 0
        self.electric_dark_correction_enabled = ''
        self.nonlinearity_correction_enabled = ''
        self.boxcar_width = 0
        self.xaxis_mode = ''
        self.number_of_pixels = 0
        self.wavelengths = []
        self.counts = []
    
    def load_file(self, filename):
        re = __import__('re')
        self.filename = filename
        with open(filename, 'r') as f:
            data = f.readlines()
            found_data = 0
            for line in data:
                line = line.rstrip()
                if found_data > 0:
                    self.add_data(line)
                if line.lower().startswith('data'):
                    self.name = line.split()[2]
		    txt = self.name 
                    re1='.*?'	# Non-greedy match on filler
		    re2='\\d+'	# Uninteresting: int
		    re3='.*?'	# Non-greedy match on filler
		    re4='\\d+'	# Uninteresting: int
		    re5='.*?'	# Non-greedy match on filler
		    re6='(\\d+)'	# Integer Number 1
		    re7='(-)'	# Any Single Character 1
		    re8='(\\d+)'	# Integer Number 2
		    re9='(-)'	# Any Single Character 2
		    re10='(\\d+)'	# Integer Number 3
		    re11='(-)'	# Any Single Character 3
		    re12='(\\d+)'	# Integer Number 4

		    rg = re.compile(re1+re2+re3+re4+re5+re6+re7+re8+re9+re10+re11+re12,re.IGNORECASE|re.DOTALL)
		    m = rg.search(txt)
		    if m:
			int1=float(m.group(1))*3600.*1000.
			c1=m.group(2)
			int2=float(m.group(3))*60.*1000.
			c2=m.group(4)
			int3=float(m.group(5))*1000.
			c3=m.group(6)
			int4=float(m.group(7))
		    self.time = int1+int2+int3+int4
                if line.lower().startswith('user'):
                    self.user = line.split()[1]
                if line.lower().startswith('spectrometer'):
                    self.spectrometer = line.split()[1]
                if line.lower().startswith('trigger'):
                    self.trigger_mode = int(line.split()[2])
                if line.lower().startswith('integration'):
                    self.integration_time = float(line.split()[3])
                if line.lower().startswith('scans'):
                    self.scans_to_average = int(line.split()[3])
                if line.lower().startswith('electric'):
                    self.electric_dark_correction_enabled = line.split()[4]
                if line.lower().startswith('nonlinearity'):
                    self.nonlinearity_correction_enabled = line.split()[3]
                if line.lower().startswith('boxcar'):
                    self.boxcar_width = int(line.split()[2])
                if line.lower().startswith('xaxis'):
                    self.xaxis_mode = line.split()[2]
                if line.lower().startswith('number'):
                    self.number_of_pixels = int(line.split()[5])
                if 'begin spectral data' in line.lower():
                    found_data = 1
        

     
    def add_data(self,line):
        """ Adding data to the wavelengths and counts lists from a given line"""
        try:
            self.wavelengths.append(float(line.split()[0]))
            self.counts.append(float(line.split()[1]))
        except ValueError:
            print 'Theres an error in the file!'
            print self.filename
            print line


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
