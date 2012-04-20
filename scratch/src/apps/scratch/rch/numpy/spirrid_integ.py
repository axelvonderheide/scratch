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
# Created on Apr 9, 2010 by: rch

import pylab as p

from stats.pdistrib.pdistrib import \
    IPDistrib, PDistrib

from numpy import \
    ogrid, frompyfunc, array, \
    linspace, hstack, vstack, vectorize, \
    apply_over_axes, sum, where, ones_like, \
    zeros_like
    
from numpy.random import random_integers 

from scipy.integrate import \
    trapz, simps, romb

from math import pi

import time

# Quantities for the response function
# and randomization
# 
E_mod = 70 * 1e+9 # Pa
sig_u = 1.25 * 1e+9 # Pa
D     = 26 * 1.0e-6 # m
A     = (D/2.0)**2 * pi
EA    = E_mod * A
xi_u  = sig_u / E_mod

def run():

    print 'DISTR XI'
    
    pxi = PDistrib( distr_choice = 'weibull_min' )
    pxi.distr_type.set( scale = 0.02, shape = 8 )
    print 'DISTR THETA'
    
    ptheta = PDistrib( distr_choice = 'uniform' )
    ptheta.distr_type.set( location = 0.0, scale = 0.01 )

    print 'DISTR LAMBDA'
    
    plambda = PDistrib( distr_choice = 'uniform' )
    plambda.distr_type.set( location = 0.0, scale = 0.2 )


    print 'PREPARING THE STATISTICAL DOMAIN'
    np = 55
    n_eps =  70
    plot_realizations = True
    
    #-------------------------------------------------------------------
    # DISCRETIZATION OF THE RANDOM VARIABLES
    #-------------------------------------------------------------------
    # this should be done in the separate object
    # called RV 
    #
    rxi     = pxi.range
    rtheta  = ptheta.range
    rlambda = plambda.range
    
    theta_ogrid = ogrid[ rxi[0]:rxi[1]:complex( 0, np ), 
                         rtheta[0]:rtheta[1]:complex( 0, np ),
                         rlambda[0]:rlambda[1]:complex( 0, np ) ]
    
    #-------------------------------------------------------------------
    # PDFS FOR INDIVIDUAL RANDOM VARIABLE DISCRETIZATIONS
    #-------------------------------------------------------------------
    # get the pdf values for xi and theta values
    # NOTE: the xi and theta have different orientations
    #       while pdf_xi goes along the 0-axis i.e. [:,0]
    #       the pdf_theta goes along the 1-axis i.e. [0,1]
    #
    pdf_xi     = pxi.distr_type.distr.pdf( theta_ogrid[0] )
    pdf_theta  = ptheta.distr_type.distr.pdf( theta_ogrid[1] )
    pdf_lambda = plambda.distr_type.distr.pdf( theta_ogrid[2] )

    #-------------------------------------------------------------------
    # PDF VALUES OVER THE DISCRETIZED RANDOM DOMAIN
    #-------------------------------------------------------------------
    # multiply the discretized pdfs to get the values of pdf for
    # all combinations of xi and theta and lambda
    # NOTE: the result of multiplying [10,1] * [1,10] array is
    #       [10,10] array
    #
    pdf_grid = pdf_xi * pdf_theta * pdf_lambda

    #------------------------------------------------------------------------------------
    # FULLY EXPANDED ORTHOGONAL GRID OF VARIABLES (including epsilon)
    #------------------------------------------------------------------------------------
    # Construct an orthogonal grid of arrays for 
    # epsilon: in the first dimension (note that empty axes 1- and 2- were
    #          using the construct eps[:,None,None] to prepare it for
    #          broadcasting (implicit copy of values for any combination of
    #          [xi,theta]
    # xi     : in the second dimension (add the first empty dimension)
    #          for broadcasting (the same xi is taken for any value of eps
    # theta  : the same as xi - therefore the list loop can be used.
    # lambda : the same as xi - therefore the list loop can be used.
    #
    eps_arr = linspace(0., 0.05, n_eps )
    
    eps = eps_arr[:,None,None,None]
    eps_theta_ogrid = [ eps ] + [ theta[None,...] for theta in theta_ogrid ]

    grid_size = reduce( lambda x,y: x * y, [ s.size for s in eps_theta_ogrid ] )
    print 'GRID SIZE: %d:' % grid_size
    if grid_size > 20000000:
        print 'That\'s too much! You have to implement it better'
        print 'or just don\'t take that many points ;-)'
        print 'you might double this when running only one version'
        return

    print '---------------------------------------------------------------'
    print 'STARTING VERSION 1 VECTORIZED RESPONSE FUNCTION EVALUATION ...'
    #---------------------------------------------------------------------------
    # VERSION 1 - SCALAR RESPONSE FUNCTION WITH SUBSEQUENT VECTORIZATION
    #---------------------------------------------------------------------------
    # evaluate the response function for all combinations of [eps,xi,theta] 
    # (see Part 1, Stochastic Modeling of Multifilament Yarns, 06)
    def q1_fn( epsilon, xi, theta, lambda_ ):
        # eps is global and eps_ is local strain
        eps_ = ( ( epsilon - theta * ( 1 + lambda_ ) ) / 
                ( ( 1 + theta ) * ( 1 + lambda_ ) ) )
        if eps_ < 0 or eps_ > xi :
            return 0.0
        else:
            return EA * eps_

    # vectorize the response function to work on arrays
    # it gets three arrays as input (strain, xi, theta, lambda%)
    #
    q1_vct_fn = vectorize( q1_fn, [float] )
        
    #----------------------------------------------------------------------------
    # CALCULATION
    #----------------------------------------------------------------------------

    start     = time.clock() # take the start time
    q1_grid    = q1_vct_fn( *eps_theta_ogrid ) # calculate (* make tuple/list to args)
    end       = time.clock() # take the end time

    print 'Duration of response function evaluation:', end - start, 's'

    if plot_realizations:
        s1 = q1_grid.shape[0]
        s2 = q1_grid.shape[1] * q1_grid.shape[2] * q1_grid.shape[3]
        q1_grid_ = q1_grid.reshape( ( s1, s2 ) )

    # multiply the response function with the pdf grid 
    # q[eps,xi,theta,lambda] * pdf[xi,theta,lambda]
    # NOTE: pdf has no eps dimension - thus, it will be broadcasted
    # to all values of eps, i.e. the same pdf_grid will be used for
    # all strains
    #
    q1_x_grid = q1_grid * pdf_grid

    # perform the integration over theta and xi and lambda_
    # NOTE: after the inner call to trapz, the resulting gets
    # reduced by one dimension, therefore, integration starts from
    # the index 3 and goes to 1
    #
    start = time.clock()
    mu_q1_arr = trapz(
                    trapz( 
                        trapz( 
                               q1_x_grid, 
                               theta_ogrid[2].flatten(), axis = 3 
                               ), 
                        theta_ogrid[1].flatten(), axis = 2 
                        ),
                    theta_ogrid[0].flatten(), axis = 1 
                )
    end = time.clock()
    
    print 'Duration of the integration:             ', end - start, 's'
    
    p.plot( eps_arr, mu_q1_arr, color = 'red', linewidth = 5 )
    
    print '---------------------------------------------------------------'
    print 'STARTING VERSION 2 RESPONSE FUNCTION USING ARRAY EXPRESSIONS ...'    
    #--------------------------------------------------------------------------
    # VERSION 2: RESPONSE FUNCTION USING ARRAY EXPRESSIONS
    #--------------------------------------------------------------------------
    #
    # implement the response function with arrays as variables
    # first extract the variable discretizations from the orthogonal grid
    #
    start     = time.clock() # take the start time
    
    eps, xi, theta, lambda_ = eps_theta_ogrid

    # NOTE: as each variable is an array oriented in different direction
    # the algebraic expressions (-+*/) perform broadcasting,. i.e. performing
    # the operation for all combinations of values. Thus, the resulgin eps
    # is contains the value of local strain for any combination of 
    # global strain, xi, theta and lambda 
    #
    eps_ = ( eps - theta * ( 1 + lambda_ ) ) / ( ( 1 + theta ) * ( 1 + lambda_ ) )

    # broadcast eps also in the xi - dimension 
    # (by multiplying with array containg ones with the same shape as xi )
    #
    eps_grid = eps_ * ones_like( xi )  
    q2_grid = EA * eps_grid
    q2_grid[ eps_grid <   0 ] = 0
    q2_grid[ eps_grid >= xi ] = 0

    end = time.clock()
    print 'Duration of response function evaluation:', end - start, 's'
    
    if plot_realizations:

        # plot all the realizations stored in sig_grid
        # first reshape the sig_grid such that it has
        # epsilons in the first index and the index s1 of the 
        # random parameter combination in the second index s2.
        # (s1 and s2 was defined above)
        #
        q2_grid_ = q2_grid.reshape( ( s1, s2 ) )
        
        # take randomly some 2000 realizations
        if s2 > 4000:
            r = random_integers( 0, high = s2, size = 4000 )
        else:
            r = range( s2 )
        for s in r:
            q2_arr = q2_grid_[:,s]
            p.plot( eps_arr, q2_arr, color = 'grey' )
        
    q2_x_grid = q2_grid * pdf_grid

    start = time.clock()    
    mu_q2_arr = trapz(
                    trapz( 
                        trapz( 
                               q2_x_grid, 
                               theta_ogrid[2].flatten(), axis = 3 
                               ), 
                        theta_ogrid[1].flatten(), axis = 2 
                        ),
                    theta_ogrid[0].flatten(), axis = 1 
                )
    end = time.clock()
    
    print 'Duration of the integration:             ', end - start, 's'
    print '--------------------------------------------------------------'

    p.plot( eps_arr, mu_q2_arr, color = 'blue', linewidth = 2 )

    p.show()

if __name__ == '__main__':

    use_profiling = False
    
    if use_profiling:
        import cProfile
        cProfile.run('run()', 'spirrid_integ_tprof' )
        
        import pstats
        p = pstats.Stats('spirrid_integ_tprof')
        p.strip_dirs()
        print 'cumulative'
        p.sort_stats('cumulative').print_stats(50)
        print 'time'
        p.sort_stats('time').print_stats(50)
    else:
        run()