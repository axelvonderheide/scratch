'''
Created on 26.06.2011

@author: axel
'''
from math import pi, e
from numpy import linspace, round , sqrt, sign , trapz
from numpy.random import rand
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt
from numpy.random import randn
from scipy.stats import norm
from enthought.traits.api import \
    Instance, Enum, Bool, on_trait_change, Int, Event, Array, Tuple, List, \
    Float, HasTraits, Float, Property, Button

from enthought.traits.api import \
    Float, Str, implements

from enthought.traits.ui.ui_traits import Image

from enthought.traits.ui.menu import OKButton, CancelButton

from enthought.traits.ui.api import \
    View, Item

from math import e, pi
from numpy import sqrt, linspace, sign, abs, cos
#from stats.spirrid.i_rf import IRF
#from stats.spirrid.rf import RF

from matplotlib import pyplot as plt

def H( x ):
    return sign( sign( x ) + 1. )

class CBShortFiber_Dev( HasTraits ):

    xi = Float( .1 )

    E_f = Float( 200e+3 )

    D_f = Float( 0.3 )

    le = Float( 8.5 )

    L_f = Float( 17.0 )

    tau = Float( 1.76 )
    f = Float( 0.03 )

    phi = Float( 0.0 )

    l = Float( 0.0 )

    theta = Float( 0.0 )

    w = Array( ctrl_range = ( 0, 0.016, 200 ) )

    def __call__( self, w, tau, L_f, D_f, E_f, le, phi, f, l, theta, xi ):
            l = l * ( 1 + theta )
            #w = w - theta * l
            T = tau * pi * D_f
            E = E_f
            A = D_f ** 2 / 4. * pi
            # debonding stage
            q_deb = -l * T + sqrt( ( l * T ) ** 2 + E * A * T * w )
            #print w
            q_deb_dev = E * A * T * 0.5 * ( ( l * T ) ** 2 + E * A * T * w ) ** ( -0.5 )
            #print q_deb_dev
            # displacement at which debonding is finished
            w0 = le * T * ( le + 2 * l ) / E_f / A

            # pulling out stage - the fiber is pulled out from the
            # side with the shorter embedded length only
            q_pull = le * T * ( ( w0 - w ) / ( ( le + 1e-15 ) - w0 ) + 1 )
            q_pull_dev = -( ( le * T ) / ( ( le + 1e-15 ) - w0 ) + 1 ) * w ** 0

            #print q_pull
            #print q_pull.shape
            q = ( q_deb * H( le * T - q_deb ) + q_pull * H( q_deb - le * T ) ) * e ** ( f * phi )
            q_dev = ( q_deb_dev * H( le * T - q_deb ) + q_pull_dev * H( q_deb - le * T ) ) * H( A * E_f * xi - q ) * e ** ( f * phi )

            # include inclination influence
            #q = q * H( q ) * e ** ( f * phi )
            #q_dev = q_dev * H( q ) * e ** ( f * phi )

            # include breaking strain
            #q = q_dev * H( A * E_f * xi - q )
            #print q
            return q_dev

if __name__ == '__main__':
    po = CBShortFiber_Dev()
    w = linspace( 0.0, 0.2, 100 )
    P = po( w, 1.76, 17.0, 0.3, 200e3, 8.5, 0.0, 0.03, 1.0, 0.02, 0.0179 )
    plt.plot( w, P )
    plt.show()



''' IDEA''
  random_field_fibers = Property( Array, depends_on = 'Length', 'height', 'width', 'Fiber volume fraction' )

    def _get_random_field_fibers( self ):
        rf = GaussRandomField( xgrid = self.x_arr , lacor = 1.7 , nsim = 1 , mean = 636 , stdev = 25.0 )
        return rf.random_field
'''
