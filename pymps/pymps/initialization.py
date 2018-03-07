
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








