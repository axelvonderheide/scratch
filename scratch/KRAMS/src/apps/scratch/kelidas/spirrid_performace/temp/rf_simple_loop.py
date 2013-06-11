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
    HasTraits, Float, Str, implements

from stats.spirrid.i_rf import \
    IRF

from stats.spirrid.rf import \
    RF

import os

from numpy import sign, linspace, array

from matplotlib import pyplot as plt

from scipy.weave import inline, converters

from types import ListType

def Heaviside( x ):
    return ( sign( x ) + 1.0 ) / 2.0

class Filament( RF ):

    implements( IRF )

    title = Str( 'brittle filament' )

    x = Float( 0.017857, auto_set=False, enter_set=True,
                distr=['weibull_min', 'uniform', 'norm'] )
    y = Float( 0.01, auto_set=False, enter_set=True, distr=['uniform', 'norm'] )
    z = Float( 0.01, auto_set=False, enter_set=True, distr=['uniform', 'norm'] )
    a = Float( 0.01, auto_set=False, enter_set=True, distr=['uniform', 'norm'] )
    C_code = '''
                for(int i=0; i < 1000000; i++)
                  q = i*x+y+z+a;
        '''

    def __call__( self, eps, x, y, z, a ):
        '''
        Implements the response function with arrays as variables.
        first extract the variable discretizations from the orthogonal grid.
        '''
        q_grid = 0
        for i in range( 0, 1000000 ):
            q_grid = i * x + y + z + a

        return q_grid

if __name__ == '__main__':
    f = Filament()

    print 'keys', f.param_keys
    print 'values', f.param_list

    print 'uniform', f.traits( distr=lambda x: x != None and 'uniform' in x )

    X = linspace( 0, 0.05, 100 )
    Y = []
    for eps in X:
        Y.append( f( eps, .017, .01, .2, 5.30929158457e-10, 70.e9 ) )
    plt.plot( X, Y, linewidth=2, color='navy' )
    plt.show()

