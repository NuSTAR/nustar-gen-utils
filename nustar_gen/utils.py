import os
import warnings
import numpy as np


def energy_to_chan(keV):
    '''
    Convert keV to NuSTAR PI channels
    PI = (keV - 1.6) / 0.04
    
    chan = energy_to_chan(3.0)
    
    
     Example
    -------
    >>> from nustar_gen.utils import energy_to_chan
    >>> chan = energy_to_chan(10.)
    >>> np.isclose(chan, 210)
    True

    '''
    return int((keV - 1.6) / 0.04)

