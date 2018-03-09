
from math import floor, ceil, sqrt

import pympscore

def estim_nband(dh: float, re:float) -> int:
    return int(ceil(re/dh))
####


def count_nzero_grid(wfun,
                     dh: float, re: float,
                     ndim = 0) -> float:
    """Count n_0 on simple grid."""
    nzero = 0.0

    nband = estim_nband(dh, re)

    if ndim <= 0:
        ndim = pympscore.GetSpaceDim()

    if ndim == 2:
        for i in range(-nband, nband+1):
            for j in range(-nband, nband+1):
                if i!=0 or j!=0:
                    x = dh * i
                    y = dh * j
                    r = sqrt(x**2 + y**2)
                    w = wfun(r, re)
                    nzero += w
    elif ndim == 3:
        for i in range(-nband, nband+1):
            for j in range(-nband, nband+1):
                for k in range(-nband, nband+1):
                    if i!=0 or j!=0 or k!=0:
                        x = dh * i
                        y = dh * j
                        z = dh * k
                        r = sqrt(x**2 + y**2 + z**2)
                        w = wfun(r, re)
                        nzero += w
    else:
        raise ValueError("Invalid ndim=%s" % ndim)

    return nzero
####

def count_lambda_grid(wfun, dh:float, re:float, ndim=0) -> float:
    """Count lambda for Laplacian"""
    lam = 0.0
    wgt = 0.0

    nband = estim_nband(dh, re)

    if ndim <= 0:
        ndim = pympscore.GetSpaceDim()

    if ndim == 2:
        for i in range(-nband, nband+1):
            for j in range(-nband, nband+1):
                if i!=0 or j!=0:
                    x = dh * i
                    y = dh * j
                    r = sqrt(x**2 + y**2)
                    w = wfun(r, re)
                    lam += w * r**2
                    wgt += w
    elif ndim == 3:
        for i in range(-nband, nband+1):
            for j in range(-nband, nband+1):
                for k in range(-nband, nband+1):
                    if i!=0 or j!=0 or k!=0:
                        x = dh * i
                        y = dh * j
                        z = dh * k
                        r = sqrt(x**2 + y**2 + z**2)
                        w = wfun(r, re)
                        lam += w * r**2
                        wgt += w
    else:
        raise ValueError("Invalid ndim=%s" % ndim)

    lam /= wgt

    return lam
####


