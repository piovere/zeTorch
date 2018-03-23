""" Contains class for reading in one run datafile.
"""

class Run(object):
    """
    Class for one run datafile's contents.
    
    This class contains every piece of information for a single data file.

    Attributes
    ----------
    filename : str
    name : str
        Name within the file
    user : str
        Username
    spectrometer : str
        Spectrometer name
    trigger_mode : int
        Triggermode of the spectrometer
    integration_time : float
    scans_to_average -- Integer of the scans to average value
    electric_dark_correction_enabled -- String
    nonlinearity_correction_enabled -- String
    boxcar_width -- Integer of the width of the boxcar
    xaxis_mode -- String
    number_of_pixels -- Integer
    wavelengths -- List of floats of the wavelengths in units of nm
    counts -- List of floats of the counts 
    
    Methods
    -------
    add_data(self,line) -- Definition to add wavelength and count data to their respective
    lists
    load_file(self, filename) -- Loads the contents of the file and assigns it to the
    respective attribute 
    """   


    
    def __init__(self):
        """Docstring
        """
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
    
    def load_file(self, filename):
        """Docstring
        """
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
        """ Adding data to the wavelengths and counts lists from a given line
        """
        try:
            self.wavelengths.append(float(line.split()[0]))
            self.counts.append(float(line.split()[1]))
        except ValueError:
            print 'Theres an error in the file!'
            print self.filename
            print line
