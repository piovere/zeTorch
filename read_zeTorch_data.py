#!/bin/python

###################################################################################################
### Name: Read zeTorch Data
### Developer: Tucker McClanahan
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
import sys
import os
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def main():
### Getting the alltabs file
### Run some checks
    directory_path, starting_wavelength, ending_wavelength = run_checks()
    
    ### Parse the files in the directory
    all_data = parse_input(directory_path)
    
    plot_data(all_data, directory_path)
    
    print 'Total number of data files: %d' % np.size(all_data)
        
    refined_data = refine_data(all_data, starting_wavelength, ending_wavelength)
    
    print 'Total number of Refined data files: %d' % np.size(refined_data)
    
    plot_refined_data(refined_data, directory_path)
    plot_data(refined_data, directory_path)

def plot_data(all_data, path):
    path = path.replace('/','_')
    pp = PdfPages(path+'_All_Data.pdf')
    
    times = range(np.size(all_data))
    lam = np.asarray(all_data[0].wavelengths)
     
    T, L = np.meshgrid(times, lam)

    data = np.zeros([np.size(lam), np.size(times)])

    for t in times:
        for l in range(np.size(lam)):
            data[l,t] = all_data[t].counts[l]
    plt.figure()
    plt.title('Plot of All Data')
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Time (AU)')
    plt.pcolormesh(L, T, data, rasterized=True, cmap='hot')
    plt.colorbar(label='Counts')
    plt.savefig(pp, format='pdf')

    pp.close()


def plot_refined_data(refined_data, path):
    path = path.replace('/','_')
    pp = PdfPages(path+'.pdf')
    
    peak_heights = []
    peak_areas = []
    peak_fwhms = []

    for run in refined_data:
        peak_heights.append(run.max)
        peak_areas.append(run.area)
        peak_fwhms.append(run.fwhm)
        
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
    plt.figure()
    plt.title('Time vs. Peak FWHM')
    plt.xlabel('Time (au)')
    plt.ylabel('Peak FWHM (counts)')
    plt.plot(range(np.size(peak_fwhms)),peak_fwhms, 'bo')
    plt.savefig(pp, format='pdf')
   
    pp.close() 

def refine_data(all_data, p_s, p_e):

    data = []
    indexes = [i for i in range(np.size(all_data[0].wavelengths)) if all_data[0].wavelengths[i]>p_s and all_data[0].wavelengths[i]<p_e]
    for run in all_data:
        ## Find the peak around 655 +- 10
        peak_lam = run.wavelengths[indexes[0]:indexes[-1]]
        peak_data = run.counts[indexes[0]:indexes[-1]]
        peak_height = np.amax(peak_data)
        if peak_height < 2000:
            continue
        
        run.max = peak_height
        run.area = np.trapz(peak_data,peak_lam)
        run.fwhm = FWHM(peak_lam, peak_data)
	
	data.append(run)

    return data


def FWHM(x,y):
    half_max = max(y)/2.
    d = np.sign(half_max - np.array(y[0:-1])) - np.sign(half_max - np.array(y[1:]))
    left_idx = np.argwhere(d > 0)[0,0]
    right_idx = np.argwhere(d<0)[0,-1]
    return x[right_idx]-x[left_idx] 


        


    
     
    
def parse_input(dirpath):
    if not dirpath.endswith('/'):
        dirpath = dirpath+'/'
    filenames = [f for f in os.listdir(dirpath) if os.path.isfile(os.path.join(dirpath, f))]
    
    data = []

    for filename in filenames:
        if filename.startswith('.'):
            continue
        with open(dirpath+filename, 'r') as f:
            erics_data = f.readlines()
            found_data = 0
            run = Run()
            run.filename = filename
            for line in erics_data:
                line = line.rstrip()
                if found_data > 0:
                    run.add_data(line)
                if line.lower().startswith('data'):
                    run.name = line.split()[2]
                if line.lower().startswith('user'):
                    run.user = line.split()[1]
                if line.lower().startswith('spectrometer'):
                    run.spectrometer = line.split()[1]
                if line.lower().startswith('trigger'):
                    run.trigger_mode = int(line.split()[2])
                if line.lower().startswith('integration'):
                    run.integration_time = float(line.split()[3])
                if line.lower().startswith('scans'):
                    run.scans_to_average = int(line.split()[3])
                if line.lower().startswith('electric'):
                    run.electric_dark_correction_enabled = line.split()[4]
                if line.lower().startswith('nonlinearity'):
                    run.nonlinearity_correction_enabled = line.split()[3]
                if line.lower().startswith('boxcar'):
                    run.boxcar_width = int(line.split()[2])
                if line.lower().startswith('xaxis'):
                    run.xaxis_mode = line.split()[2]
                if line.lower().startswith('number'):
                    run.number_of_pixels = int(line.split()[5])
                if 'begin spectral data' in line.lower():
                    found_data = 1
                
        data.append(run)

    return data


class Run(object):
    
    def __init__(self):
        self.filename = ''
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
        self.max_655 = 0.0
        self.area_655 = 0.0
        self.fwhm_655 = 0.0
         
    def add_data(self,line):
        try:
            self.wavelengths.append(float(line.split()[0]))
            self.counts.append(float(line.split()[1]))
        except ValueError:
            print 'Theres an error in the file!'
            print self.filename
            print line


def run_checks():
    # Run checks on alltabs_file_input
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
