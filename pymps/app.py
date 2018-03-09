
import sys
from typing import List

import numpy as NP
from numpy import linalg as LA

import pymps
from pymps.number_density import count_nzero_grid, count_lambda_grid
from pymps.visualization import write_fluiddata_vtu
from pymps.initialization import gen_simple_block, load_prof_file
from pymps.arrayutil import *

import pympscore


class ConstWeight(pympscore.StandardWeightFunction):
    '''Constant weight function. Always return 1 if within influence.
    NOTE this extends a C-class'''
    def __init__(self):
        return super().__init__()
    def eval(self, r, re):
        return 1.0 if r<re else 0.0
    def __call__(self, r, re):
        return self.eval(r, re);
####

class App2d(object):
    '''2D fluid application'''

    def __init__(self):
        # bootstrap
        self.ndim = 2
        pympscore.SetSpaceDim(ndim)
        print('spacedim=%d' % pympscore.GetSpaceDim())

        # physical
        init_phys(self)

        # numerical
        init_comp(self)

    ####

    def init_phys(self):
        '''Init physical'''
        self.phys_rho = 1.0e3
        self.phys_cs = 10.0
        self.phys_pref = self.phys_rho * self.phys_cs**2
        self.phys_nu = 1.0e-6
        self.phys_grav = NP.array([0.0, -9.81], dtype=float)
        return
    ####

    def init_comp(self):
        '''Init computational'''

        # discretization
        dh = 0.01
        self.dh = dh
        self.re_dens = dh * 2.1
        self.re_grad = dh * 2.1
        self.re_visc = dh * 4.1
        self.re_pres = dh * 4.1
        self.re_max = max(re_dens, re_grad, re_visc, re_pres)

        # weight function
        wfun = pympscore.StandardWeightFunction()
        wgradfun = wfun
        self.wfun = wfun
        self.wgradfun = wgradfun

        #
        nzero = count_nzero_grid(wfun, dh, re_dens)
        nzero_part = count_nzero_grid((lambda r,re: 1 if r<re else 0), dh, re_dens)
        nzero_grad = count_nzero_grid(wgradfun, dh, re_grad)
        nzero_visc = count_nzero_grid(wfun, dh, re_visc)
        nzero_pres = count_nzero_grid(wfun, dh, re_pres)
        nlam_visc = count_lambda_grid(wfun, dh, re_visc)
        nlam_pres = count_lambda_grid(wfun, dh, re_pres)
        nsurf = nzero * 0.97
        relax_pres = 0.2
        nlam_pres *= relax_pres
        print('nzero=%g' % nzero)
        print('nzero_part=%g' % nzero_part)
        print('nzero_grad=%g' % nzero_grad)
        print('nzero_visc=%g' % nzero_visc)
        print('nzero_pres=%g' % nzero_pres)
        print('lambda_visc=%g' % nlam_visc)
        print('lambda_pres=%g' % nlam_pres)

        self.nzero = nzero
        self.nzero_part = nzero_part
        self.nzero_grad = nzero_grad
        self.nzero_visc = nzero_visc
        self.nzero_pres = nzero_pres
        #
        self.nlam_visc = nlam_visc
        self.nlam_pres = nlam_pres
        #
        self.nsurf = nsurf

        return
    ####

    def init_data(self):

        return
    ####


    def main_wc(args: List[str]):
        '''Weakly compressible algo.'''

        return
    ####


####


def main2d():
    '''2D run'''
    
    # bootstrap
    ndim = 2
    pympscore.SetSpaceDim(ndim)
    print('spacedim=%d' % pympscore.GetSpaceDim())

    # physical
    rho = 1.0e3
    cs = 10.0
    pref = rho * cs**2
    nu = 1.0e-1
    grav = NP.array([0.0, -9.81])

    # discretization
    dh = 0.01
    re_dens = dh * 2.1
    re_grad = dh * 2.1
    re_visc = dh * 4.1
    re_pres = dh * 4.1
    re_max = max(re_dens, re_grad, re_visc, re_pres)

    # weight function
    wfun = pympscore.StandardWeightFunction()
    wgradfun = wfun
    print(wfun)

    nzero = count_nzero_grid(wfun, dh, re_dens)
    nzero_part = count_nzero_grid(ConstWeight(), dh, re_dens)
    nzero_grad = count_nzero_grid(wgradfun, dh, re_grad)
    nzero_visc = count_nzero_grid(wfun, dh, re_visc)
    nzero_pres = count_nzero_grid(wfun, dh, re_pres)
    nlam_visc = count_lambda_grid(wfun, dh, re_visc)
    nlam_pres = count_lambda_grid(wfun, dh, re_pres)
    nsurf = nzero * 0.97
    relax_pres = 0.2
    nlam_pres *= relax_pres
    print('nzero=%g' % nzero)
    print('nzero_part=%g' % nzero_part)
    print('nzero_grad=%g' % nzero_grad)
    print('nzero_visc=%g' % nzero_visc)
    print('nzero_pres=%g' % nzero_pres)
    print('lambda_visc=%g' % nlam_visc)
    print('lambda_pres=%g' % nlam_pres)


    # data 
    prob_type = -1
    prob_type = 0
    fluid_data = pympscore.FluidData(ndim)
    npart = 0
    if prob_type == 0:
        init_pos = []
        #init_pos += gen_simple_block(dh, [100,100], [-0.5,-0.5])
        #init_pos += gen_simple_block(dh, [100,50], [-1.0,-0.25])
        #init_pos += gen_simple_block(dh, [50,100], [0.0,-0.5])
        init_pos += gen_simple_block(dh, [100, 50], [-0.5, 0.0])
        init_pos += gen_simple_block(dh, [200, 10], [-1.0, -0.1])
        npart = len(init_pos)

        # 
        fluid_data.numParticle = npart
        fluid_data.alloc()
        NP.copyto(fluid_data.position(), init_pos)

        # mask walls
        msk: NP.ndarray = (fluid_data.position()[:,1] <= 0)
    else:
        prof = load_prof_file('input.grid', ndim)
        npart = prof['npart']

        #
        fluid_data.numParticle = npart
        fluid_data.alloc()
        NP.copyto(fluid_data.position(), prof['pos'])
        NP.copyto(fluid_data.surfflag(), prof['type'])

        msk: NP.ndarray = (fluid_data.surfflag() != 0)

    # domain
    domain = pympscore.Domain()
    domain.vlo = [-2.5, -2.5]
    domain.vhi = [2.5, 2.5]
    print(domain)

    # cache
    re_big = re_max * 1.05
    bucket = pympscore.BucketCache(domain)
    bucket.define(re_big)
    print(bucket)
    bucket.cachePositions(fluid_data.position())

    # neighbor
    neigh_table = pympscore.NeighborTable()
    neigh_table.cutoff = re_big
    neigh_table.resize(fluid_data.numParticle)
    neigh_table.buildFromCache(bucket, fluid_data.position())
    print("table size=%d" % len(neigh_table))

    # expose data access
    numdens: NP.ndarray = fluid_data.numdens()
    surfflag: NP.ndarray = fluid_data.surfflag()
    pympscore.calcNumDens(wfun, re_dens, neigh_table, numdens)
    print("dens_min=%f" % numdens.min())
    print("dens_max=%f" % numdens.max())

    #
    coords: NP.ndarray = fluid_data.position()
    vel: NP.ndarray = fluid_data.velocity()
    pres: NP.ndarray = fluid_data.pressure()
    # assign some values
    vel.fill(0)
    pres.fill(0)

    # local buffers
    grad = newVectorArray(npart, ndim)
    visc = newVectorArray(npart, ndim)

    #
    iplot = 0
    write_fluiddata_vtu("./output/homu%04d"%iplot, fluid_data)
    iplot += 1

    #sys.exit(0)

    time = 0.0
    dt = 0.001
    max_step = 1000
    for step in range(1, max_step+1):

        # build neighbor list
        bucket.cachePositions(coords)
        neigh_table.buildFromCache(bucket, coords)

        # compute number density
        pympscore.calcNumDens(wfun, re_dens, neigh_table, numdens)
        
        # tag surface part
        pympscore.tagSurfFlag(nsurf, numdens, surfflag)

        # compute viscosity
        for dir in range(ndim):
            coef = 2.0*ndim/(nlam_visc*nzero_visc)
            pympscore.calcLaplace(wfun, re_visc, coef, neigh_table, vel[:,dir], visc[:,dir])

        # compute pressure
        if False:
            # explicit use EOS
            pympscore.calcPresEos(numdens, pres, pref, nzero)

            pympscore.calcSymGrad(wgradfun, re, nzero_grad, neigh_table, coords, pres, grad)

            # update
            acc = (-1.0/rho) * grad + nu * visc + grav
            #acc[:,0] += grav[0]
            #acc[:,1] += grav[1]
            vel += acc * dt
            vel[msk,:] = 0
            coords += vel * dt
        else:
            # semi-implicit 

            # trial move
            acc = nu * visc + grav
            vel += acc * dt
            vel[msk,:] = 0
            coords += vel * dt

            neigh_table.updateInfo(coords)

            #
            pympscore.calcNumDens(wfun, re_dens, neigh_table, numdens)
            pympscore.tagSurfFlag(nsurf, numdens, surfflag)

            #
            rhs = (1.0/dt) / nzero * (numdens-nzero)

            #
            coef = 2.0*ndim / (nlam_pres * nzero_pres)
            coef *= (dt/rho)
            pympscore.solvePres(wfun, re_pres, 
                                coef, 
                                neigh_table, surfflag, 
                                rhs, pres)
            
            # clip pressure
            pres[pres<0] = 0

            # 
            pympscore.calcSymGrad(wfun, re_grad, nzero_grad, neigh_table, coords, pres, grad)

            # update
            acc = (-1.0/rho) * grad
            acc[msk,:] = 0
            vel += acc * dt
            coords += acc * dt**2

        # advance time
        time += dt

        if step%10 == 0:
            print('step=%d, time=%f' % (step,time))
            write_fluiddata_vtu("./output/homu%04d"%iplot, fluid_data)
            iplot += 1

    return 
####

if __name__ == "__main__":

    main2d()




####




