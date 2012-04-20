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
# Created on Jun 14, 2010 by: rch

from enthought.traits.api import \
    HasTraits, Float, Str, implements, Range, Property, cached_property, Array

from math import pi

from numpy import \
    sign, sqrt, linspace

from quaducom.resp_func.cb_clamped_fiber import CBClampedFiber

from stats.spirrid.i_rf import \
    IRF

from stats.spirrid.rf import \
    RF

from matplotlib import pyplot as plt

def Heaviside( x ):
    return ( sign( x ) + 1.0 ) / 2.0

class CBClampedX( CBClampedFiber ):

    implements( IRF )

    # @todo - define them as a range
    # - define embedded loops for compiled eps_loop
    # - visualization of nd response? - mayavi, cutting, slicing in spirrid_view?
    # 

    x = Property( Array, depends_on = 'L1, L2, l', ctrl_range = ( 0.0, 200.0, 100 ) )
    @cached_property
    def _get_x( self ):
        return linspace( 0, self.L1 + self.L2 + self.l, 100 )

    w = Float( auto_set = False, enter_set = True,
               ctrl_range = ( 0.0, 0.05, 100 ) )


    def __call__( self, w, tau, l, D_f, E_f, theta, xi, phi, L1, L2 ):
        '''Calculate the stress transfer 
        length associated with the current crack width.
        '''
        T = tau * phi * D_f * pi
        q = super( CBClampedX, self ).__call__( w, tau, l, D_f, E_f, theta, xi, phi, L1, L2 )
        x_l = self.x - l
        return ( q - ( T * x_l * Heaviside( x_l ) ) ) * Heaviside( q - T * x_l )

if __name__ == '__main__':

    cbx = CBClampedX()
    P = cbx( 0.5, 2.5, 0.0, 26e-3, 72e3, 0.01, 1, 1, 1, 0.5 )
    plt.plot( cbx.x, P )
    plt.show()
