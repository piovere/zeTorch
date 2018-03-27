""" Contains class for reading in one run datafile.
"""
import re 
import numpy as np
from scipy.interpolate import interp1d

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
    
    def calibrate_data(self, calibration_file, lamp_det_file):
        cal_data = np.loadtxt(calibration_file, delimiter=',', skiprows = 1)
        lamp_data = np.loadtxt(lamp_det_file, skiprows= 14)
        lam = np.asarray(self.wavelengths)
        counts = np.asarray(self.counts)
        
        min_lam = max(cal_data[:,0].min(), lamp_data[:,0].min(), np.min(lam))
        max_lam = min(np.max(cal_data[:,0]), np.max(lamp_data[:,0]), np.max(lam))
       
        indexes = [i for i in range(np.size(cal_data[:,0])) if cal_data[i,0] > min_lam and cal_data[i,0]< max_lam] 
        slashed_cal_data = cal_data[indexes,:]
        indexes = [i for i in range(np.size(lamp_data[:,0])) if lamp_data[i,0] > min_lam and lamp_data[i,0]< max_lam] 
        slashed_lamp_data = lamp_data[indexes,:]
        indexes = [i for i in range(np.size(lam)) if lam[i] > min_lam and lam[i]< max_lam] 
        slashed_lam = lam[indexes]
        slashed_counts = counts[indexes]

        print(slashed_cal_data.shape)

        cal_interp = interp1d(slashed_cal_data[:,0], slashed_cal_data[:,1], bounds_error=False, fill_value=1.0)
        lamp_interp = interp1d(slashed_lamp_data[:,0], slashed_lamp_data[:,1], bounds_error=False, fill_value=1.0)
        meas_interp = interp1d(slashed_lam, slashed_counts, bounds_error=False, fill_value=1.0)
        
        print(slashed_lam[-1])
        print(slashed_cal_data[-1,0])
        print(slashed_lamp_data[-1,0])

        self.calibrated_data = np.dstack((slashed_lam, cal_interp(slashed_lam) * slashed_counts / lamp_interp(slashed_lam)))[0]
        np.savetxt('might_be_calibrated.txt', self.calibrated_data)
        
        
    def load_file(self, filename):
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
