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
# Created on Sep 8, 2011 by: rch

from enthought.traits.api import implements, Str
from stats.spirrid import \
    SPIRRID, RV, RF, IRF
from stats.spirrid import SPIRRIDLAB
import numpy as np
import numexpr as ne
import math
from scipy.interpolate import interp1d

def Heaviside(x):
        ''' Heaviside function '''
        return x >= 0#( np.sign( x ) + 1.0 ) / 2.0

def Heaviside_ne(x):
    ''' Heaviside function '''
    return ne.evaluate("x >= 0")#( np.sign( x ) + 1.0 ) / 2.0


if __name__ == '__main__':

    #===========================================================================
    # Response function
    #===========================================================================
    class fiber_tt_5p_np(RF):
        ''' Response function of a single fiber '''
        implements(IRF)

        title = Str('brittle filament')

        def __call__(self, eps, lambd, xi, E_mod, theta, A):
            '''
            Implements the response function with arrays as variables.
            first extract the variable discretizations from the orthogonal grid.
            '''
            # NOTE: as each variable is an array oriented in different direction
            # the algebraic expressions (-+*/) perform broadcasting,. i.e. performing
            # the operation for all combinations of values. Thus, the resulgin eps
            # is contains the value of local strain for any combination of 
            # global strain, xi, theta and lambda 
            #

            eps_ = (eps - theta * (1 + lambd)) / ((1 + theta) * (1 + lambd))

            # cut off all the negative strains due to delayed activation
            # 
            eps_ *= Heaviside(eps_)

            # broadcast eps also in the xi - dimension 
            # (by multiplying with array containing ones with the same shape as xi )
            #
            eps_grid = eps_ * Heaviside(xi - eps_)

            # cut off all the realizations with strain greater than the critical one.
            # 
            # eps_grid[ eps_grid >= xi ] = 0

            # transform it to the force
            # 
            q_grid = E_mod * A * eps_grid

            return q_grid

        C_code = '''
                double eps_ = ( eps - theta * ( 1 + lambd ) ) /
                                 ( ( 1 + theta ) * ( 1 + lambd ) );
                // Computation of the q( ... ) function
                if ( eps_ < 0 || eps_ > xi ){
                    q = 0.0;
                }else{
                      q = E_mod * A * eps_;
                }
            '''

    class fiber_tt_5p_ne(RF):
        ''' Response function of a single fiber '''
        implements(IRF)

        title = Str('brittle filament')

        def __call__(self, eps, lambd, xi, E_mod, theta, A):
            '''
            Implements the response function with arrays as variables.
            first extract the variable discretizations from the orthogonal grid.
            '''
            # NOTE: as each variable is an array oriented in different direction
            # the algebraic expressions (-+*/) perform broadcasting,. i.e. performing
            # the operation for all combinations of values. Thus, the resulgin eps
            # is contains the value of local strain for any combination of 
            # global strain, xi, theta and lambda 
            #

#            eps_ = ne.evaluate("(eps - theta * (1 + lambd)) / ((1 + theta) * (1 + lambd))")
#
#            # cut off all the negative strains due to delayed activation
#            # 
#            eps_ *= Heaviside(eps_)
#
#            # broadcast eps also in the xi - dimension 
#            # (by multiplying with array containing ones with the same shape as xi )
#            #
#            tmp = Heaviside(ne.evaluate("xi - eps_"))
#            eps_grid = ne.evaluate("eps_ * tmp")
#
#            # cut off all the realizations with strain greater than the critical one.
#            # 
#            # eps_grid[ eps_grid >= xi ] = 0
#
#            # transform it to the force
#            # 
#            q_grid = ne.evaluate("E_mod * A * eps_grid")

            return ne.evaluate("E_mod * A *((eps - theta * (1 + lambd)) / ((1 + theta) * (1 + lambd))) * (((eps - theta * (1 + lambd)) / ((1 + theta) * (1 + lambd)))>=0) *  ((xi-((eps - theta * (1 + lambd)) / ((1 + theta) * (1 + lambd))))>=0)")

        C_code = '''
                double eps_ = ( eps - theta * ( 1 + lambd ) ) /
                                 ( ( 1 + theta ) * ( 1 + lambd ) );
                // Computation of the q( ... ) function
                if ( eps_ < 0 || eps_ > xi ){
                    q = 0.0;
                }else{
                      q = E_mod * A * eps_;
                }
            '''

    D = 26 * 1.0e-6 # m
    A = (D / 2.0) ** 2 * math.pi

    # set the mean and standard deviation of the two random variables
    la_mean, la_stdev = 0.0, 0.2
    xi_mean, xi_stdev = 0.019027, 0.0022891
    E_mean, E_stdev = 70.0e+9, 15.0e+9
    th_mean, th_stdev = 0.0, 0.01
    A_mean, A_stdev = A * 0.3, 0.7 * A

    do = 'norm'

    if do == 'general':

        # set the mean and standard deviation of the two random variables
        la_mean, la_stdev = 0.0, 0.2
        xi_mean, xi_stdev = 0.019027, 0.0022891
        E_mean, E_stdev = 70.0e+9, 15.0e+9
        th_mean, th_stdev = 0.0, 0.01
        A_mean, A_stdev = A * 0.3, 0.7 * A

        # construct the normal distributions and get the methods
        # for the evaluation of the probability density functions
        g_la = RV('uniform', la_mean, la_stdev)
        g_xi = RV('norm', xi_mean, xi_stdev)
        g_E = RV('uniform', E_mean, E_stdev)
        g_th = RV('uniform', th_mean, th_stdev)
        g_A = RV('uniform', A_mean, A_stdev)

        mu_ex_file = 'fiber_tt_5p_30.txt'
        delimiter = ','

    elif do == 'uniform':

        # set the mean and standard deviation of the two random variables
        la_mean, la_stdev = 0.0, 0.2
        xi_mean, xi_stdev = 0.01, 0.02
        E_mean, E_stdev = 70.0e+9, 15.0e+9
        th_mean, th_stdev = 0.0, 0.01
        A_mean, A_stdev = A * 0.3, 0.7 * A

        # construct the uniform distributions and get the methods
        # for the evaluation of the probability density functions
        g_la = RV('uniform', la_mean, la_stdev)
        g_xi = RV('uniform', xi_mean, xi_stdev)
        g_E = RV('uniform', E_mean, E_stdev)
        g_th = RV('uniform', th_mean, th_stdev)
        g_A = RV('uniform', A_mean, A_stdev)

        mu_ex_file = 'fiber_tt_5p_40_unif.txt'
        delimiter = ' '

    elif do == 'norm':

        # set the mean and standard deviation of the two random variables
        la_mean, la_stdev = 0.1, 0.02
        xi_mean, xi_stdev = 0.019027, 0.0022891
        E_mean, E_stdev = 70.0e+9, 15.0e+9
        th_mean, th_stdev = 0.005, 0.001
        A_mean, A_stdev = 5.3e-10, 1.0e-11

        # construct the normal distributions and get the methods
        # for the evaluation of the probability density functions
        g_la = RV('norm', la_mean, la_stdev)
        g_xi = RV('norm', xi_mean, xi_stdev)
        g_E = RV('norm', E_mean, E_stdev)
        g_th = RV('norm', th_mean, th_stdev)
        g_A = RV('norm', A_mean, A_stdev)

        mu_ex_file = 'fiber_tt_5p_40_norm.txt'
        delimiter = ' '

    # discretize the control variable (x-axis)
    e_arr = np.linspace(0, 0.04, 1000)

    #===========================================================================
    # Randomization
    #===========================================================================


    s_ne = SPIRRID(q = fiber_tt_5p_ne(),
                e_arr = e_arr,
                n_int = 10,
                tvars = dict(lambd = g_la, xi = g_xi, E_mod = g_E, theta = g_th, A = g_A),
                )
    s_np = SPIRRID(q = fiber_tt_5p_np(),
                e_arr = e_arr,
                n_int = 10,
                tvars = dict(lambd = g_la, xi = g_xi, E_mod = g_E, theta = g_th, A = g_A),
                )

    print 'np', s_np.exec_time
    print 'ne', s_ne.exec_time
