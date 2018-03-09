
import numpy as NP

from math import floor, ceil


def gen_simple_block(dh: float, num, vlo):

    ndim = len(num)

    data = []

    if ndim == 2:
        nx,ny = num
        for j in range(ny):
            for i in range(nx):
                x = dh * (i+0.5) + vlo[0]
                y = dh * (j+0.5) + vlo[1]
                data.append((x,y))
    elif ndim == 3:
        nx,ny,nz = num
        for k in range(nz):
            for j in range(ny):
                for i in range(nx):
                    x = dh * (i+0.5) + vlo[0]
                    y = dh * (j+0.5) + vlo[1]
                    z = dh * (k+0.5) + vlo[2]
                    data.append((x,y,z))
    else:
        raise ValueError("Invalid dimension")

    return data
####

def load_prof_file(filename: str, ndim=3) -> dict:
    '''Load old PROF file'''
    prof = {}
    prof['type'] = []
    prof['pos'] = []
    with open(filename) as fp:
        #skip time
        fp.readline()

        npart = int(fp.readline())
        
        prof['npart'] = npart

        for i in range(npart):
            line:str = fp.readline()
            fields = line.strip().split()
            prof['type'].append(int(fields[0]))
            prof['pos'].append([float(x) for x in fields[1:1+ndim]])

    return prof
####






