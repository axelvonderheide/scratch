'''
Created on Jun 17, 2010

@author: rostislav
'''

from enthought.traits.api import HasTraits, Array, Float
from scipy.special import gamma
from scipy.optimize import leastsq
from numpy import linspace, array, max, sqrt, min
from matplotlib import pyplot as plt
from scipy.stats import weibull_min

class ParamsExtraction( HasTraits ):

    # tested parameters    
    tested_l = Array
    tested_fu = Array
    l0 = Float( 20. ) # mm, reference length
    
    # define function for calculating the Weibull scaling
    def powerlaw( self, length, scale, shape ):
        return scale * gamma( 1. + 1. / shape ) * ( self.l0 / length ) ** ( 1. / shape )
    
    def error_func( self, p ):
        return self.tested_fu - self.powerlaw( self.tested_l, p[0], p[1] )
    
    def weibull_params( self ):
        iparams = [110., 3.0]
        return leastsq( self.error_func, iparams, full_output = 1 )
    
    E = Float( 70.e9 )
    l = Float( 200. )
    # stress-strain curve of a perfect clamped bundle consisting of
    # fibers with two parameters Weibull strength distribution given
    # by the evaluated parameters
    def values( self, eps ):
        out = self.weibull_params()
        scale0 = out[0][0]
        shape = out[0][1]
        scale = ( scale0 * ( self.l0 / self.l ) ** ( 1. / shape ) ) / self.E * 1.e6
        print shape, scale
        return eps * self.E * ( 1. - weibull_min.cdf( eps, shape, scale = scale ) )

if __name__ == '__main__':

#################################
######### PLOTTING ##############
#################################
    
    pe = ParamsExtraction( tested_l = array( [11., 20., 50., 130.] ), \
                          tested_fu = array( [2.949, 2.427, 2.227, 2.027] ) )
    
    def weibull_params_fit():
        out = pe.weibull_params()
        params = out[0]
        covar = out[1]
        scale = params[0]
        shape = params[1]
        shapeErr = sqrt( covar[0][0] )
        scaleErr = sqrt( covar[1][1] )
    
        x = linspace( 9., 135., 100 )
        y = pe.powerlaw( x, scale, shape )
        
        plt.subplot( 2, 1, 1 )
        plt.plot( pe.tested_l, pe.tested_fu, 'ro' )
        plt.plot( x, y )
        shape_ytick = min( y ) + 0.9 * ( max( y ) - min( y ) )
        scale_ytick = min( y ) + 0.8 * ( max( y ) - min( y ) )
        xtick = min( x ) + 0.5 * max( x ) - min( x )
        plt.text( xtick, shape_ytick, 'shape: %5.2f +/- %5.2f' % ( shape, shapeErr ) )
        plt.text( xtick, scale_ytick, 'scale for $l_0$ = ' + str( pe.l0 ) + 'mm: %5.2f +/- %5.2f' % ( scale, scaleErr ) )
        plt.title( 'best fit Weibull scaling' )
        plt.xlabel( 'length' )
        plt.ylabel( 'strength' )
        plt.xlim( x[0] - ( x[-1] - x[0] ) * 0.05, x[-1] + ( x[-1] - x[0] ) * 0.05 )
        plt.ylim( y[-1] + ( y[-1] - y[0] ) * 0.05, y[0] - ( y[-1] - y[0] ) * 0.05 )
        
        plt.subplot( 2, 1, 2 )
        plt.loglog( x, y )
        plt.loglog( pe.tested_l, pe.tested_fu, 'ro' )
        plt.xlabel( 'log length' )
        plt.ylabel( 'log strength' )
        plt.xlim( x[0] - ( x[-1] - x[0] ) * 0.01, x[-1] + ( x[-1] - x[0] ) * 0.05 )
        plt.ylim( y[-1] + ( y[-1] - y[0] ) * 0.05, y[0] - ( y[-1] - y[0] ) * 0.05 )
        plt.show()
        
    def bundle():
        eps = linspace( 0, 0.03, 100 )
        sigma = pe.values( eps )
        plt.plot( eps, sigma, linewidth = 2 )
        plt.title( 'stress-strain of an asymptotic bundle' )
        plt.xlabel( 'strain [-]' )
        plt.ylabel( 'stress [MPa]' )
        plt.show()
    
    #bundle()
    weibull_params_fit()
    
    
    
        
