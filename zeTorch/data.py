"""Contains class for data files. Maybe should eventually sublass Spectrum?
"""

import numpy as np


class Data(object):
    """DOCSTRING
    """
    def __init__(self, file=None):
        """DOCSTRING
        """
        self._file = file
    
    def load(self, file=None):
        if file is None and self._file is None:
            raise FileNotFoundError("You must provide a filename")
        
        if file is None:
            file = self._file

        self._data = np.loadtxt(file, skiprows=14)

    @property
    def data(self):
        return self._data
    
    @property
    def wavelengths(self):
        return self._data[:,0]
    
    @property
    def intensity(self, )
