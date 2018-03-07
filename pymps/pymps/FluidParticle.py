
import numpy

#from . import getSpaceDim

class FluidParticle(object):
    """description of class"""

    # fields annotation
    parttype: numpy.ndarray
    position: numpy.ndarray
    velocity: numpy.ndarray
    pressure: numpy.ndarray

    #
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        ndim = 2
        print(ndim)

        self.parttype = None
        self.position = None
        self.velocity = None
        self.pressure = None

    ####

#### FluidParticle

class Timer(object):
    '''Description'''

    # fields annotation


    pass
#### Timer





