'''
Created on May 11, 2010

@author: rostislav
'''

''' module for plotting the strain, displacement and stress in a yarn held by two clamps '''

from enthought.traits.api import HasTraits, Float, Property, Tuple, Array, Enum, Range
from enthought.traits.ui.api import View, Item
from mathkit.mfn.mfn_line.mfn_line import MFnLineArray
from matplotlib import pyplot as plt
from numpy import linspace, hstack, sqrt, tanh, sinh, cosh, array, frompyfunc, \
    trapz, sign
from scipy.optimize import brentq

def H( x ):
    return sign( sign( x ) + 1 )

class ClampConst( HasTraits ):

    # distance between the Einspannklemme and Kraftschlussklemme
    l1 = Float( 0.137, modified = True, auto_set = False, enter_set = True )

    #tested length
    lt = Float( 0.1, modified = True, auto_set = False, enter_set = True )

    # length Kraftschlussklemme
    Lk = Float( 0.05, modified = True, auto_set = False, enter_set = True )

    # yarn module of elasticity
    E = Float( 72e9, modified = True, auto_set = False, enter_set = True )

    # yarn cross-sectional area
    A = Float( 0.89e-6, modified = True, auto_set = False, enter_set = True )

    # current tensile force
    Ft = Range( 0., 1000., modified = True, auto_set = False, enter_set = True )

    # compression dependent shear flow per length - Kraftschlussklemme 
    qk = Float( 2200., modified = True, auto_set = False, enter_set = True )

    # pressure Kraftschlussklemme
    Pk = Float( 6.0, modified = True, auto_set = False, enter_set = True )

    # additional force on the Kraftschlussklemme
    Fa = Range( 10.0, 700.0, modified = True, auto_set = False, enter_set = True )

    # prestressing force
    P0 = Float( 100., modified = True, auto_set = False, enter_set = True )

    values = Tuple( Array, Array )

    def strains( self ):
        mfn = MFnLineArray()
        mfn.xdata, mfn.ydata = self.values
        strains_fn = frompyfunc( mfn.get_diff, 1, 1 )
        strains = strains_fn( mfn.xdata )
        strains[0] = strains[1]
        strains[-2] = strains[-1]
        return strains

    def get_values( self ):
        l1 = self.l1
        lt = self.lt
        A = self.A
        E = self.E
        Ft = self.Ft
        Fa = self.Fa
        Pk = self.Pk
        qk = self.qk
        Lk = self.Lk
        P0 = self.P0

        # before the add-clamp moves axially 
        xdata = linspace( 0, l1 + Lk + lt, 200 )
        def ydata1( x ):
            y0 = min((l1 + Lk - x)*qk*Pk,P0,Lk*qk*Pk)
            y1 = min((Ft-Fa),(Ft-Fa)-qk*Pk*(x-l1))
            return max(y0,y1)

        py_y1 = frompyfunc( ydata1, 1, 1 )
        y1 = py_y1( xdata )

        def ydata2( x ):
            y2 = Ft
            y3 = (x - Lk - l1)*qk*Pk + Ft
            return min(y2,y3)

        py_y2 = frompyfunc( ydata2, 1, 1 )
        y2 = py_y2( xdata )

        y1 = y1 * H( y1 )
        y2 = y2 * H( y2 )

        f = frompyfunc( lambda x, y: max( x, y ) , 2, 1 )
        y = f( y1, y2 )

        self.values = xdata, y

        return self.values



    traits_view = View( Item( 'l1', label = 'fixation to clamping distance [m]' ),
                       Item( 'lt', label = 'tested length [m]' ),
                       Item( 'Lk', label = 'clamp length [m]' ),
                       Item( 'E', label = 'yarn elasticity modulus [N/m2]' ),
                       Item( 'A', label = 'yarn cross-section [m2]' ),
                       Item( 'Ft', label = 'current tensile force [N]' ),
                       Item( 'qk', label = 'frictional shear force per pressure [m-1]' ),
                       Item( 'Pk', label = 'applied force to the clamp [N]' ),
                       Item( 'Fa', label = 'additional axial force [N]' ),
                       Item( 'P0', label = 'prestressing force [N]' ),
                       )

if __name__ == '__main__':
    cs = ClampConst()
    x, y = cs.get_values()
    plt.plot( x, y )
    plt.show()
