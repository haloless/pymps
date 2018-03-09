
import pympscore

import numpy as NP
import pyevtk

def write_fluiddata_vtu(prefix: str, fd: pympscore.FluidData):
    ntotal = fd.numParticle
    ndim = fd.numDim

    zeros = NP.zeros(ntotal) if ndim<3 else None

    pos = fd.position()
    xpos = pos[:,0]
    ypos = pos[:,1]
    zpos = pos[:,2] if ndim==3 else zeros

    data = {}
    data['pressure'] = fd.pressure()
    data['numdens'] = fd.numdens()
    data['velocity'] = (fd.velocity()[:,0], fd.velocity()[:,1], 
                        fd.velocity()[:,2] if ndim==3 else zeros)
    data['surfflag'] = fd.surfflag()
    #data['type'] = fd.

    pyevtk.hl.pointsToVTK(prefix, xpos,ypos,zpos, data)

    return
####






