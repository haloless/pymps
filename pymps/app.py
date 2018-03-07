
import sys

import numpy as NP
from numpy import linalg as LA

import pymps
import pymps.number_density
from pymps.visualization import write_fluiddata_vtu
from pymps.initialization import gen_simple_block

import pympscore







def main2d():
    # bootstrap
    ndim = 2
    pympscore.SetSpaceDim(ndim)
    print('spacedim=%d' % pympscore.GetSpaceDim())

    #
    rho = 1.0e3
    cs = 10.0
    pref = rho * cs**2

    # 
    dh = 0.01
    re_scale = 3.1
    re = dh * re_scale

    # weight function
    wfun = pympscore.ModifiedWeightFunction()
    wgradfun = pympscore.ModifiedGradientWeightFunction()

    nzero = pymps.number_density.count_nzero_grid(wfun, dh, re)
    nzero_grad = pymps.number_density.count_nzero_grid(wgradfun, dh, re)
    print('nzero=%f' % nzero)
    print('nzero_grad=%f' % nzero_grad)

    # domain
    domain = pympscore.Domain()
    domain.vlo = [-2.5,-2.5]
    domain.vhi = [2.5,2.5]
    print(domain)

    # 
    init_pos = []
    #init_pos += gen_simple_block(dh, [100,100], [-0.5,-0.5])
    init_pos += gen_simple_block(dh, [100,50], [-1.0,-0.25])
    init_pos += gen_simple_block(dh, [50,100], [0.0,-0.5])
    npart = len(init_pos)

    fluid_data = pympscore.FluidData()
    fluid_data.numDim = ndim
    fluid_data.numParticle = npart
    
    # 
    fluid_data.alloc()
    NP.copyto(fluid_data.position(), init_pos)

    #
    bucket = pympscore.BucketCache(domain)
    bucket.define(re)
    print(bucket)
    bucket.cachePositions(fluid_data.position())

    #
    neigh_table = pympscore.NeighborTable()
    neigh_table.cutoff = re
    neigh_table.resize(fluid_data.numParticle)
    neigh_table.buildFromCache(bucket, fluid_data.position())
    print("table size=%d" % len(neigh_table))

    neigh_list = neigh_table[0] # access by index
    for neigh in neigh_list.asList():
        print("j=%d, rij=%f" % (neigh.index, neigh.distance))

    #
    numdens: NP.ndarray = fluid_data.numdens()
    pympscore.calcNumDens(wfun, neigh_table, numdens)

    print("dens_min=%f" % numdens.min())
    print("dens_max=%f" % numdens.max())


    #
    coords: NP.ndarray = fluid_data.position()
    vel: NP.ndarray = fluid_data.velocity()
    pres: NP.ndarray = fluid_data.pressure()

    # assign some values
    vel.fill(0)
    vel[coords[:,0]<=0,0] = 0.5
    vel[coords[:,0]>0,0] = -0.5
    #coords[coords[:,0]<=0,1] += 0.2
    #coords[coords[:,0]>0,1] += -0.2

    pres.fill(0)

    grad = NP.empty((npart,ndim), dtype=NP.float64, order='F')

    iplot = 0
    write_fluiddata_vtu("./output/homu%04d"%iplot, fluid_data)
    iplot += 1

    time = 0.0
    dt = 0.0002
    max_step = 1000
    for step in range(1, max_step+1):

        # 
        bucket.cachePositions(coords)
        neigh_table.buildFromCache(bucket, coords)

        #
        pympscore.calcNumDens(wfun, neigh_table, numdens)

        pympscore.calcPresEos(numdens, pres, pref, nzero)

        pympscore.calcSymGrad(wgradfun, re, nzero_grad, neigh_table, coords, pres, grad)


        # update
        vel += (-dt/rho) * grad 
        coords += vel * dt


        time += dt
        print('step=%d, time=%f' % (step,time))

        if step%10 == 0:
            write_fluiddata_vtu("./output/homu%04d"%iplot, fluid_data)
            iplot += 1

    return 
####

if __name__ == "__main__":

    main2d()




####




