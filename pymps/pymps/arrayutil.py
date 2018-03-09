"""Helper to create data arrays"""

import numpy as NP

import pympscore

def newScalarArray(numpart: int, data_type=float) -> NP.ndarray:
    '''Create new scalar array'''
    return NP.empty(numpart, dtype=data_type)
####

def newVectorArray(numpart:int, ndim:int = 0, data_type=float, data_order:str='F') -> NP.ndarray:
    """Create new vector array. 
    If NDIM not given, use MPS space dimension.
    The array is by default Fortan-style, to be consistent with Eigen.
    """
    if ndim <= 0:
        ndim = pympscore.GetSpaceDim()
    return NP.empty((numpart,ndim), dtype=data_type, order=data_order)
####


