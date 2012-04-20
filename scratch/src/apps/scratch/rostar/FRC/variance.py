'''
Created on Jan 4, 2011

@author: rostislav
'''

from math import pi
from matplotlib import pyplot as plt
from numpy import arccos, sin, linspace, mean, var, sum, trapz, arange
from numpy.random import rand
from scipy.stats import norm, weibull_min

class Test():

    N_shape = 4.5
    N_scale = 10
    N_loc = 10
    X_shape = 7
    X_scale = 12
    X_loc = 5
    sims = 500
    
    def N( self, x ):
        return weibull_min.pdf( x, self.N_shape, scale = self.N_scale, loc = self.N_loc )
    
    def N_ppf( self, n ):
        return weibull_min.ppf( n, self.N_shape, scale = self.N_scale, loc = self.N_loc )
    
    def X( self, x ):
        return weibull_min.pdf( x, self.X_shape, scale = self.X_scale, loc = self.X_loc )

    def X_ppf( self, n ):
        return weibull_min.ppf( n, self.X_shape, scale = self.X_scale, loc = self.X_loc )
    
if __name__ == '__main__':
    
    x = linspace( 0, 40, 1000 )
    
    var_sim = []
    var_eval = []
    no_sim = linspace( 100, 10000, 30 )
    no_sim = [10000000]
    for points in no_sim:
        t = Test()
        t.sims = points
        X = t.X( x )
        N = t.N( x )
        sim = t.sims
        result = []
        filaments = t.N_ppf( rand( sim ) )
        for i, s in enumerate( filaments ):
            resp = sum( t.X_ppf( rand( s - 1 ) ) )
            result.append( resp )
        mean_N = weibull_min( t.N_shape, scale = t.N_scale, loc = t.N_loc ).stats( 'm' )
        var_N = weibull_min( t.N_shape, scale = t.N_scale, loc = t.N_loc ).stats( 'v' )
        mean_X = weibull_min( t.X_shape, scale = t.X_scale, loc = t.X_loc ).stats( 'm' )
        var_X = weibull_min( t.X_shape, scale = t.X_scale, loc = t.X_loc ).stats( 'v' )
        var_D = mean_N * var_X + mean_X ** 2 * var_N
        var_sim.append( var( result ) )
        var_eval.append( var_D )
    
    print var_sim
    print var_eval
    plt.plot( no_sim, var_sim )
    plt.plot( no_sim, var_eval )
    #print 'computed variance = ', var_D  
    #print 'variance for', sim, 'simulated data =', var( result )
    #plt.hist( result, normed = True, bins = 50 )
    plt.show()
    
