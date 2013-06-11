#-------------------------------------------------------------------------------
#
# Copyright (c) 2009, IMB, RWTH Aachen.
# All rights reserved.
#
# This software is provided without warranty under the terms of the BSD
# license included in simvisage/LICENSE.txt and may be redistributed only
# under the conditions described in the aforementioned license.  The license
# is also available online at http://www.simvisage.com/licenses/BSD.txt
#
# Thanks for using Simvisage open source!
#
# Created on Mar 23, 2010 by: rch

from ibvpy.mats.mats_explore import \
    MATSExplore

from ibvpy.mats.mats2D.mats2D_explore import MATS2DExplore
    
from ibvpy.mats.mats2D.mats2D_cmdm.mats2D_cmdm import \
    MATS2DMicroplaneDamage
    
from ibvpy.mats.matsXD.matsXD_cmdm.matsXD_cmdm_phi_fn import \
     PhiFnStrainHardeningLinear, PhiFnStrainSoftening

import math

import numpy

from enthought.mayavi import mlab
 
def run_trial_simulation():

    phi_fn = PhiFnStrainHardeningLinear( alpha = 0.4, beta = 0.4, Elimit = 1.0 )
    #phi_fn = PhiFnStrainSoftening()

    explorer = MATSExplore( dim = MATS2DExplore( mats_eval = \
                                                 MATS2DMicroplaneDamage(n_mp = 30, 
                                                                        phi_fn = phi_fn) ))

    explorer.tloop.bcond_list[0].max_strain = 0.006
    explorer.tloop.bcond_list[0].alpha_rad  = 0.0 # math.pi / 4.0
    
    u = explorer.tloop.eval()

    tsig = explorer.tloop.tstepper.rtrace_mngr[ 'time - sig_norm' ]
    tsig.refresh()
    tsig.configure_traits()
    
    print tsig.trace.xdata
    print tsig.trace.ydata
#    print 'u', u
    
    return

@mlab.show
def scan_strain_space( max_strain = 0.0002 ):

    rho = 0.03
    E_f = 70000
    E_m = 24000  
    E_c = E_m * ( 1 - rho ) + E_f * rho


    R = max_strain
    n_R = 30
    in_R = 30j
    
    Alpha_min =  - math.pi / 2.0 * 0.5
    Alpha_max =    math.pi / 2.0 * 1.5
    n_Alpha = 19
    in_Alpha = 19j
    
    phi_fn = PhiFnStrainHardeningLinear(  E_m = E_m, E_f = E_f, rho = rho,
                                          sigma_0 = 5.0, 
                                          alpha = 0.4, beta = 0.4, Elimit = 1.0 )
    #phi_fn = PhiFnStrainSoftening()

    explorer = MATSExplore( dim = MATS2DExplore( mats_eval = \
                                                 MATS2DMicroplaneDamage(E = E_c, nu = 0.25,
                                                                        n_mp = 30, 
                                                                        phi_fn = phi_fn) ))

    explorer.tloop.bcond_list[0].max_strain = R
    explorer.tloop.tline.step = 1.0 / (n_R-1)
    teps = explorer.tloop.tstepper.rtrace_mngr[ 'strain - strain' ]
    teps.idx_y = 3
#    r, alpha = numpy.mgrid[0:R:in_R, Alpha_min:Alpha_max: in_Alpha ]
#    print 'r', r
#    print 'alpha', alpha
#    
#    x = r * numpy.cos( alpha )
#    y = r * numpy.sin( alpha )

    def f( alpha_rad ):
        explorer.tloop.reset()        
        print 'calculating', alpha_rad
        explorer.tloop.bcond_list[0].alpha_rad = alpha_rad
        u = explorer.tloop.eval()
        tsig = explorer.tloop.tstepper.rtrace_mngr[ 'time - sig_norm' ]
        #teps = explorer.tloop.tstepper.rtrace_mngr[ 'strain - strain' ]        
        tsig.refresh()
        teps.refresh()
        return teps.trace.xdata, teps.trace.ydata, tsig.trace.ydata            

    alpha_arr = numpy.linspace( Alpha_min, Alpha_max, n_Alpha )
    print 'alpha_arr', alpha_arr
    
    sfunc = numpy.frompyfunc( f, 1, 3 )
    xarr_list, yarr_list, zarr_list  = sfunc( alpha_arr )
    print 'x',xarr_list
    print 'y',yarr_list
    print 'z',zarr_list
    
    x = numpy.array( numpy.vstack( xarr_list ), dtype = 'float_' ).T
    y = numpy.array( numpy.vstack( yarr_list ), dtype = 'float_' ).T
    z = numpy.array( numpy.vstack( zarr_list ), dtype = 'float_' ).T

    maxx = numpy.max( numpy.abs( x ) )
    maxy = numpy.max( numpy.abs( y ) )
    maxz = numpy.max( z ) / 0.5
    
    x /= maxx
    y /= maxy
    z /= maxz
    
    print 'x', x
    print 'y', y
    print 'z', z
    
    s = mlab.mesh(x, y, z, colormap="bone" )
    return

if __name__ == '__main__':
    scan_strain_space( 0.006 )
    #run_trial_simulation()
