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
# Created on Jun 2, 2010 by: rch


from enthought.traits.api import \
    HasTraits, Float, Str, implements, Range

from math import pi, e

from numpy import sign, linspace, array, cos, sqrt, argmax, hstack, max, zeros, ones

from matplotlib import pyplot as plt

from stats.spirrid.i_rf import \
    IRF

from stats.spirrid.rf import \
    RF

def Heaviside( x ):
    return ( sign( x ) + 1.0 ) / 2.0


class ConstantFrictionFiniteFiber( RF ):

    implements( IRF )

    title = Str( 'pull-out with constant friction' )

    fu = Float( 1200.0e6, auto_set=False, enter_set=True,
                distr=['weibull_min'] )

    qf = Float( 1500., auto_set=False, enter_set=True,
                distr=['uniform', 'norm'] )

    L = Float( 0.02, auto_set=False, enter_set=True, distr=['uniform'] )

    A = Float( 5.30929158457e-10, auto_set=False, enter_set=True, distr=['uniform', 'weibull_min'] )
    E_mod = Float( 70.0e9, auto_set=False, enter_set=True, distr=['uniform'] )
    phi = Range( 0., pi, auto_set=False, enter_set=True, distr=['cos_distr'] )
    z = Float( 0., auto_set=False, enter_set=True, distr=['uniform'] )
    f = Float( 0.01, auto_set=False, enter_set=True, distr=['uniform'] )

    C_code = '''
            double w = eps;
            double Le = L / 2. - z / cos( phi );
            double w_deb = exp( f * phi ) * qf * pow(Le,2.0) / E_mod / A;
            double P_deb_full = sqrt( w * E_mod * A * qf ) * exp( f * phi );
            double P_deb;
            
            // Heaviside
            if ( Le < 0 || P_deb_full > fu * A || w > w_deb ){
                P_deb = 0;
            }else{
                P_deb =P_deb_full;
            }
            
            double P_pull_x = ( Le * qf - Le * qf / ( Le - w_deb ) * ( w - w_deb ) ) * exp( f * phi );
            double P_pull;
            
            // Heaviside 
            if ( P_pull_x < 0 || w_deb > w ){
                P_pull = 0;
            }else{
                P_pull = P_pull_x;
            }
            
            // Computation of the q( ... ) function
            q = P_deb + P_pull;
        '''

#    def __call__( self, w, fu, qf, L, A, E_mod, z, phi, f ):
#
#        Le = -z 
#        Le = Le / cos( phi )
#        Le = Le + L / 2. 
#        w_deb = f * phi  
#        w_deb = e ** ( w_deb )
#        w_deb = w_deb * qf 
#        w_deb = w_deb * Le ** 2.0
#        w_deb = w_deb / E_mod 
#        w_deb = w_deb / A
#        P_deb_full = w * E_mod  
#        P_deb_full = P_deb_full * A 
#        P_deb_full = P_deb_full * qf 
#        P_deb_full = sqrt( P_deb_full )
#        P_deb_full = P_deb_full * e ** ( f * phi )
#        H1 = Heaviside( fu * A - P_deb_full ) 
#        P_deb_full = P_deb_full * H1
#        del H1
#        H2 = Heaviside( w_deb - w ) 
#        P_deb_full = P_deb_full * H2
#        del H2
#        H3 = Heaviside( Le )
#        P_deb_full = P_deb_full * H3
#        del H3
#        P_pull_x = Le * qf 
#        P_pull_x = P_pull_x - Le 
#        P_pull_x = P_pull_x * qf 
#        P_pull_x = P_pull_x / ( Le - w_deb ) 
#        P_pull_x = P_pull_x * ( w - w_deb )  
#        P_pull_x = P_pull_x * e ** ( f * phi )
#        P_pull_x = P_pull_x * Heaviside( P_pull_x ) 
#        P_pull_x = P_pull_x * Heaviside( w - w_deb )
#        return P_deb_full + P_pull_x

#    def __call__( self, w, fu, qf, L, A, E_mod, z, phi, f ):
#
#        Le = zeros( max( L.shape ) * max( z.shape ) * max( phi.shape ) )
#        Le.reshape( tuple( array( L.shape ) * array( z.shape ) * array( phi.shape ) ) )
#        Le -= z 
#        Le /= cos( phi )
#        Le += L / 2.
#        w_deb = e ** ( f * phi ) * qf * Le ** 2.0 / E_mod / A
#        P_deb_full = sqrt( w * E_mod * A * qf ) * e ** ( f * phi )
#        P_deb = P_deb_full * Heaviside( fu * A - P_deb_full ) * Heaviside( w_deb - w ) * Heaviside( Le )
#        P_pull_x = ( Le * qf - Le * qf / ( Le - w_deb ) * ( w - w_deb ) ) * e ** ( f * phi )
#        P_pull = P_pull_x * Heaviside( P_pull_x ) * Heaviside( w - w_deb )
#        return P_deb + P_pull

    def __call__( self, w, fu, qf, L, A, E_mod, z, phi, f ):

        Le = L / 2. - z / cos( phi )
        w_deb = e ** ( f * phi ) * qf * Le ** 2.0 / E_mod / A
        P_deb_full = sqrt( w * E_mod * A * qf ) * e ** ( f * phi )
        P_deb = P_deb_full * Heaviside( fu * A - P_deb_full ) * Heaviside( w_deb - w ) * Heaviside( Le )
        P_pull_x = ( Le * qf - Le * qf / ( Le - w_deb ) * ( w - w_deb ) ) * e ** ( f * phi )
        P_pull = P_pull_x * Heaviside( P_pull_x ) * Heaviside( w - w_deb )
        return P_deb + P_pull

class ConstantFrictionAndFreeLength( RF ):
    '''
    '''

    implements( IRF )

    title = Str( 'pull-out with constant friction and free length ' )

    tau = Float( 8, auto_set=False, enter_set=True,
                distr=['uniform'] )

    # free length
    l = Float( 1, auto_set=False, enter_set=True,
                distr=['uniform', 'norm'] )

    E = Float( 70e9, auto_set=False, enter_set=True,
               distr=['uniform'] )

    A = Float( 5.30929158457e-10, auto_set=False, enter_set=True,
               distr=['uniform', 'weibull_min'] )

    # waviness in strains
    slack = Float( 0.1, auto_set=False, enter_set=True,
               distr=['uniform'] )

    def __call__( self, u, tau, l, E, A, slack ):
        return - l * ( 1 + slack ) * tau * Heaviside( u - l * ( slack ) ) + \
                + sqrt( ( l * ( 1 + slack ) * tau ) ** 2 \
                + 2 * E * A * ( u - ( l * slack ) ) * Heaviside( u - l * ( slack ) ) )
                

if __name__ == '__main__':
    cf = ConstantFrictionFiniteFiber()
    cf_fl = ConstantFrictionAndFreeLength()
    X = linspace( 0, 0.0012, 100 )
    Y_fl = cf_fl( X, 8.0, 0.01, 210.e9, 0.004, 0.01 )
    plt.plot( X, Y_fl, linewidth=2 )
    plt.show()
