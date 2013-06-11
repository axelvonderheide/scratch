'''
Created on Oct 19, 2010

@author: rostislav

Module for statistical evaluation of single filament test.
Assumed is the two parameters Weibull distribution.
'''

from numpy import array, mean, var, linspace, log
from scipy.optimize import fmin
from matplotlib import pyplot as plt
from scipy.stats import weibull_min
from scipy.stats import norm


####### TEST DATA #######
# machine stiffness: 200 cN at strains 0,28%

eps_s = array( [ 1.78104763, 3.06000589, 2.19944303, 2.06883411] )
#eps_s = array ( [2.0169252 , 2.674148] )
#eps_s = array( [6.5011566124120588, 6.4219613493812515, 5.1165304445274558, 6.8611176748112239, 8.532373930184395, 4.9878882809128005, 9.1801420354491547, 8.9165608842634416, 4.6730088888888881, 8.375569674284213] )


####### EVALUATION ########

def histogram( data ):
    plt.hist( data, 7, normed = True, label = 'raw data', color = 'lightgrey' )

def moment_method( data ):
    mean_ = mean( data )
    var_ = var( data, ddof = 1 )
    print '#### moment method ####'
    print 'mean = ', mean_
    print 'var = ', var_
    params = fmin( moment_w, [2., 5.], args = ( mean_, var_ ) )
    print 'Weibull shape = ', params[0]
    print 'Weibull scale = ', params[1]
    # plot the Weibull fit based on the moment method
    e = linspace( 0., 0.3 * ( max( data ) - min( data ) ) + max( data ), 100 )
    #plt.plot( e, weibull_min.pdf( e, params[0], scale = params[1], loc = 0. ),
         #color = 'blue', linewidth = 2, label = 'moment method' )

def moment_w( params, *args ):
    mean = args[0]
    var = args[1]
    weibull = weibull_min( params[0], scale = params[1], loc = 0. )
    r = abs( weibull.stats( 'm' ) - mean ) + abs( weibull.stats( 'v' ) - var )
    return r

# maximum likelihood estimation for two parameter Weibull distribution
# max likelihood method l to be maximized
def maxlike( v, *args ):
    data = args
    shape, scale = v
    r = sum( log( weibull_min.pdf( data, shape, scale = scale, loc = 0. ) ) )
    return - r

# optimizing method for finding the parameters
# of the fitted Weibull distribution
def maximum_likelihood( data ):
    params = fmin( maxlike, [2., 5.], args = ( data ) )
    moments = weibull_min( params[0], scale = params[1], loc = 0. ).stats( 'mv' )
    print ' #### max likelihood #### '
    print 'mean = ', moments[0]
    print 'var = ', moments[1]
    print 'Weibull shape = ', params[0]
    print 'Weibull scale = ', params[1]
    # plot the Weibull fit according to maximum likelihood 
    e = linspace( 0., 0.3 * ( max( data ) - min( data ) ) + max( data ), 100 )
    #plt.plot( e, weibull_min.pdf( e, params[0], scale = params[1], loc = 0. ),
         #color = 'red', linewidth = 2, label = 'max likelihood' )
    plt.xlabel( 'Risspannung $\sigma_{c,u}$' , fontsize = 20 )
    plt.ylabel( 'PDF', fontsize = 20 )

if __name__ == '__main__':
    # set the data to evaluate
    data = eps_s
    #data = eps_l

    histogram( data )
    moment_method( data )
    maximum_likelihood( data )
    e = linspace( 0., 0.3 * ( max( data ) - min( data ) ) + max( data ), 100 )
    plt.plot( e, norm( 6.31, 1.47 ).pdf( e ) )
    plt.legend( loc = 'best' )
    plt.title( 'Approximation der Rissspannung mit Normalverteilung ', fontsize = 20 )
    plt.show()
