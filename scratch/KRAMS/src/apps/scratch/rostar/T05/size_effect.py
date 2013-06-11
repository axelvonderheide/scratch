'''
Created on Sep 22, 2009

@author: rostislav
'''

from enthought.traits.api import HasTraits, Float, Int, Property
from numpy import linspace, argmax
from numpy.random import rand
from math import exp
from scipy.stats import weibull_min
from matplotlib import pyplot as plt 


class SizeEffect( HasTraits ):
    '''
    Size effect depending on the yarn length
    '''
    
    l_b = Float( 0.01, auto_set = False, enter_set = True, # [m]
                  desc = 'yarn total length',
                  modified = True )
    m_f = Float( 4.5, auto_set = False, enter_set = True, # [-]
                desc = 'Weibull shape parameter for filaments',
                modified = True )
    Nf = Int( 24000, auto_set = False, enter_set = True, # [-]
                desc = 'number of filaments in yarn',
                modified = True )
    l_rho = Float( 0.02, auto_set = False, enter_set = True, # [m]
                desc = 'autocorrelation length for fiber strength dispersion',
                modified = True )
    s_rho = Float( 2500, auto_set = False, enter_set = True, # [m]
            desc = 'scale parameter for autocorrelation length',
            modified = True )

    # these parameters are called plot, but the model should not 
    # interfere with the view (@todo resolve)
    l_plot = Float( 80., auto_set = False, enter_set = True, # [m]
                desc = 'maximum yarn length',
                modified = True )

    min_plot_length = Float( 0.0001, auto_set = False, enter_set = True, # [m]
                desc = 'minimum yarn length',
                modified = True )
    
    n_points = Int( 100, auto_set = False, enter_set = True,
                desc = 'points to plot',
                modified = True )

    # autocorrelation length function
    def fl( self, l ):
        return ( self.l_rho / ( self.l_rho + l ) ) ** ( 1. / self.m_f )
    '''second option'''
    #    return (l/self.l_rho + self.l_rho/(self.l_rho + l))**(-1./self.m_f)
    
    # scale parameter depending on length
    def s( self, l ):
        return self.fl( l ) * self.s_rho
    
    def mu_f( self, l ):
        return weibull_min( self.m_f, scale = self.s( l ), loc = 0.0 ).stats( 'm' )

    def mu_b( self, l ):
        return self.m_f ** ( -1.0 / self.m_f ) * self.s( l ) * exp( -1.0 / self.m_f )

class FilBun( SizeEffect ):
    
    nf = Int( 100 )
    K = Float( 72e9 * 0.89e-6 )
    weib = Property
    def _get_weib( self ):
        return weibull_min( self.m_f, scale = self.s( 0.2 ) / self.K )
    def filaments( self ):
        no = rand( self.nf )
        strains = self.weib.ppf( no )
        return strains
        

if __name__ == '__main__':
    def se():
        se = SizeEffect()
        l = linspace( 0.0005, 50, 100000 )
        mu_f = se.mu_f( l )
        mu_b = se.mu_b( l )
        chob = se.mu_b( se.l_b )
        plt.plot( l, mu_f, linewidth = 2, color = 'grey', label = 'filament strength' )
        plt.plot( l, mu_b, linewidth = 2, color = 'black', label = 'bundle strength' )
        plt.plot( [se.l_b, 50], [chob, chob], linewidth = 2, color = 'red', label = 'yarn strength' )
        plt.plot( [se.l_b, se.l_b], [10, chob], linestyle = 'dotted',
                  linewidth = 1, color = 'black' )
        plt.xscale( 'log' )
        plt.yscale( 'log' )
        plt.legend( loc = 'best' )
        plt.xlabel( 'specimen length log(l) $[m]$', size = 'large' )
        plt.ylabel( 'strength log($\sigma_u$) $[N/mm^2]$' )
        plt.xlim( ( 0.001, max( l ) * .9 ) )
        plt.ylim( ( min( mu_b ) * 0.95, max( mu_f ) * 1.2 ) )
    def fil():
        fb = FilBun()
        strains = fb.filaments()
        x = linspace( 0., max( strains ) * 1.2, 1000 )
        for s in strains:
            plt.plot( [0., s, s], [0., fb.K * s, 0.], color = 'grey' )
        plt.plot( [0., strains[-1], strains[-1]], [0., fb.K * strains[-1], 0.],
                 color = 'grey', label = 'single filaments' )
        plt.plot( fb.weib.stats( 'm' ), fb.K * fb.weib.stats( 'm' ), 'ro' )
        plt.plot( [0., fb.weib.stats( 'm' )], 2 * [fb.K * fb.weib.stats( 'm' )],
                  color = 'red', linestyle = '--', label = 'mean filament strength' )
        bundle = fb.K * x * ( 1. - fb.weib.cdf( x ) )
        plt.plot( x, bundle, linewidth = 2, color = 'black',
                 label = 'bundle' )
        idmx = argmax( bundle )
        plt.plot( x[idmx], max( bundle ), 'bo' )
        plt.plot( [0., x[idmx]], [max( bundle ), max( bundle )],
                  linestyle = '--', color = 'blue', label = 'bundle strength' )
        plt.legend( loc = 'best' )
        plt.xlabel( 'strain [-]', size = 'large' )
        plt.ylabel( 'stress [$N/mm^2$]', size = 'large' )
        plt.title( 'l-d curves for gauge length 200 mm' )
        plt.xlim( ( 0., max( x ) * 1.1 ) )
        plt.ylim( ( 0., max( strains * fb.K ) * 1.1 ) )

    
    #fil()
    se()
    plt.show()
    
    
