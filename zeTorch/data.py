"""Contains class for data files. Maybe should eventually sublass Spectrum?
"""

import numpy as np


class Data(object):
    """DOCSTRING
    """
    def __init__(self, file=None):
        """DOCSTRING
        """
        self._filename = file
    
    def load(self, file=None):
        if file is None and self._file is None:
            raise FileNotFoundError("You must provide a filename")
