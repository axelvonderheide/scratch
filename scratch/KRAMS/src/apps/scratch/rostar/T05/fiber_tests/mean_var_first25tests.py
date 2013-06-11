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


####### TEST DATA #######
# machine stiffness: 200 cN at strains 0,28%

eps_short = array( [2.87, 3.71, 3.11, 2.62, 2.64, 3.68,
                   3.09, 3.30, 3.57, 2.07, 3.28, 2.43,
                   3.43, 2.46, 2.82, 3.53, 3.32, 3.48,
                   3.93, 2.65, 2.18, 2.55, 3.26, 3.85,
                   3.20] )

eps_long = array( [2.97, 2.28, 2.33, 2.68, 2.77, 2.72, 2.28,
                   2.93, 2.52, 2.26, 2.51, 2.42, 2.69, 2.19,
                   2.51, 2.72, 2.14, 2.8, 2.67, 2.99, 2.09,
                   2.3, 2.39, 2.88, 2.88 ] )

P_long = array( [91.07, 103.9, 78.84, 82.24, 85.72, 95.3, 109.78,
                 85.42, 95.34, 73.37, 83.72, 74.46, 100.5, 85.39,
                 87.94, 106.79, 74.52, 88.11, 100.98, 105.74, 83.51,
                 140.14, 81.67, 111.75, 92.15, ] )

P_short = array( [72.47, 128.15, 85.31, 72.41, 73.92, 93.24, 129.75,
                  89.45, 110.07, 57.96, 89.81, 67.5, 105.48, 97.43,
                  87.46, 116.69, 98.26, 87.72, 117.75, 81.37, 73.06,
                  117.37, 86.72, 110.87, 110.23 ] )

# subtracting the machine stiffness from the measured data
eps_s = eps_short - P_short / 200. * 0.28
eps_l = eps_long - P_long / 200. * 0.28

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
    plt.plot( e, weibull_min.pdf( e, params[0], scale = params[1], loc = 0. ),
         color = 'blue', linewidth = 2, label = 'moment method' )

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
    plt.plot( e, weibull_min.pdf( e, params[0], scale = params[1], loc = 0. ),
         color = 'red', linewidth = 2, label = 'max likelihood' )
    plt.xlabel( 'strain' )
    plt.ylabel( 'PDF' )

if __name__ == '__main__':
    # set the data to evaluate
    data = eps_s
    #data = eps_l

    histogram( data )
    moment_method( data )
    maximum_likelihood( data )
    plt.legend( loc = 'best' )
    plt.show()
