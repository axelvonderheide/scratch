'''
Created on 3.12.2010

Response function for a double sided pullout of a continuous filament
with different friction at both sides.

@author: Q
'''

from enthought.traits.api import \
    HasTraits, Float, Str, implements, Int

from math import pi, e

from numpy import sign, linspace, array, cos, sqrt, argmax, hstack, max

from matplotlib import pyplot as plt

from stats.spirrid.i_rf import \
    IRF

from stats.spirrid.rf import \
    RF

def H( x ):
    return ( sign( x ) + 1.0 ) / 2.0


class DoublePullout( RF ):

    implements( IRF )

    title = Str( 'double yarn pullout' )

    xi = Float( 0.014, auto_set = False, enter_set = True,
                distr = ['weibull_min', 'uniform'] )
    tau1 = Float( 6., auto_set = False, enter_set = True,
                distr = ['uniform', 'norm'] )
    tau2 = Float( 6., auto_set = False, enter_set = True,
                distr = ['uniform', 'norm'] )
    l = Float( 0.01, auto_set = False, enter_set = True,
              distr = ['uniform'] )
    A = Float( 1., auto_set = False,
              enter_set = True, distr = ['uniform', 'weibull_min'] )

    d = Float( 26e-3, auto_set = False, input = True,
          enter_set = True, distr = ['uniform', 'weibull_min'] )

    E_mod = Float( 1., auto_set = False, enter_set = True,
                  distr = ['uniform'] )
    theta = Float( 0.01, auto_set = False, enter_set = True,
                  distr = ['uniform', 'norm'] )
    # embedded length 1
    L1 = Float( 2., auto_set = False, enter_set = True,
          distr = ['uniform'] )
    # embedded length 2
    L2 = Float( 3., auto_set = False, enter_set = True,
              distr = ['uniform'] )
    # number of filaments
    Nf = Int( 1723, auto_set = False, enter_set = True,
              distr = ['uniform'] )
    # bond quality
    phi = Float( 1., auto_set = False, enter_set = True,
                  distr = ['uniform', 'norm'] )

    w = Float( auto_set = False, enter_set = True,
               ctrl_range = ( 0.0, 1.0, 10 ) )

    C_code = '''
        '''

    def __call__( self, w, tau1, tau2, l, d, A, E_mod, theta, xi, phi, L1, L2, Nf ):
        l = l * ( 1 + theta )
        w = w - theta * l
        Tau1 = tau1 * phi * d * pi
        Tau2 = tau2 * phi * d * pi
        P_ = ( -tau1 * tau2 * l + ( tau1 * tau2 * ( tau2 * tau1 * l ** 2 + 2 * w * H( w ) *
                                                    E_mod * A * tau2 + 2 * w * H( w ) * E_mod * A * tau1 ) ) ** ( 1. / 2. ) ) / ( tau2 + tau1 )
        if L1 * tau1 < L2 * tau2:
            P_ = P_ * H( Tau1 * L1 - P_ ) + Tau1 * L1 * H( P_ - Tau1 * L1 )
        else:
            P_ = P_ * H( Tau2 * L2 - P_ ) + Tau2 * L2 * H( P_ - Tau2 * L2 )
        P = P_ * H( A * E_mod * xi - P_ )
        return P * Nf

if __name__ == '__main__':
    dp = DoublePullout()
    X = linspace( 0.0, 0.25, 100 )
    Y = dp( X, 2.5, 2.5, .01, 26e-3, pi * 26e-3 ** 2 / 4., 70.0e3, 0.01, 0.014, 1., 2., 3., 1723 )
    plt.plot( X, Y, linewidth = 2 )
    plt.show()
