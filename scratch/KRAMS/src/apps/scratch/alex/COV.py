'''
Created on 16.06.2011

@author: rrypl
'''

from numpy import array, mean, var, sqrt, linspace, log, sum
from matplotlib import pyplot as plt
from scipy.stats import norm
from scipy.optimize import fmin
from math import e

F0 = array( [107.4, 107.9, 107.9, 99.0, 101.7, 96.7] )

F90 = array( [106.3, 114.4, 110.0] )

PO12c = array( [32.2, 35.7, 33.9, 33.6, 37.2, 31.8, 32.1] )

ST112c = array( [82.6, 77.3, 76.0, 80.8, 70.3, 72.8, 73.0] )

ST212c = array( [143.1, 137.0, 140.0, 143.7, 133.7, 126.3, 127.0] )

BT0 = array( [11.9, 11.8, 12.4, 12.2] )
BT0 = BT0 * 1.15 / 4


BT90 = array( [13.6, 12.6, 14.1, 12.6] )
BT90 = BT90 * 1.15 / 4
print BT90

# select experimental values
values = BT90
kn = 1.83

mu = mean( values )
stdev = sqrt( var( values ) )

print 'mean =', mu
print 'stdev =', stdev
print 'COV =', stdev / mu

Vx = stdev / mu
if Vx < 0.1:
    Vx = 0.1

my = mean( log( values ) )
sy = sqrt( log( Vx ** 2 + 1 ) )
Xk = e ** ( my - kn * sy )

print 'my =', my
print 'sy =', sy
print 'Xk =', Xk
print 'Xd =', Xk / 1.5
#print 'sigma =', Xk / 1.5 / 0.04 / 0.16 / 1000, 'MPa'

#------------ ADDITIONAL CRAP ------------------
# evaluation of the 5% quantile from distributions based entirely on the experiments

#def histogram( data ):
#    plt.hist( data, 3, normed = True, label = 'raw data', color = 'lightgrey' )
#
#def moment_norm( params, *args ):
#    mean = args[0]
#    var = args[1]
#    normal = norm( loc = params[0], scale = params[1] )
#    r = abs( normal.stats( 'm' ) - mean ) + abs( normal.stats( 'v' ) - var )
#    return r
#
#def moment_method( data ):
#    mean_ = mean( data )
#    var_ = var( data, ddof = 1 )
#    params = fmin( moment_norm, [2., 5.], args = ( mean_, var_ ) )
#    # plot the Gauss normal fit based on the moment method
#    distr = norm( loc = params[0], scale = params[1] )
#    e = linspace( mean_ - 10 * sqrt( var_ ), mean_ + 10 * sqrt( var_ ), 200 )
#    plt.plot( e, distr.pdf( e ),
#         color = 'blue', linewidth = 2, label = 'stat. moment method' )
#    # plot of the 5% quantile
#    quantile = distr.ppf( 0.05 )
#    plt.plot( 2 * [quantile], [0, max( distr.pdf( e ) )], color = 'blue', lw = 2, label = '5% qunatile =' + ' %.1f N' % quantile )
#
#def maxlike( v, *args ):
#    data = args
#    loc, scale = v
#    r = sum( log( norm.pdf( data, loc = loc, scale = scale ) ) )
#    return - r
#
## optimizing method for finding the parameters
## of the fitted normal distribution
#def maximum_likelihood( data ):
#    params = fmin( maxlike, [2., 5.], args = ( data ) )
#    # plot the Gauss normal fit according to maximum likelihood 
#    distr = norm( loc = params[0], scale = params[1] )
#    e = linspace( params[0] - 10 * sqrt( params[1] ), params[0] + 10 * sqrt( params[1] ), 200 )
#    plt.plot( e, distr.pdf( e ),
#         color = 'red', linewidth = 2, label = 'max likelihood method' )
#    quantile = distr.ppf( 0.05 )
#    plt.plot( 2 * [quantile], [0, max( distr.pdf( e ) )], color = 'red', lw = 2, label = '5% qunatile =' + ' %.1f N' % quantile )
#    plt.xlabel( 'Force [kN]' )
#    plt.ylabel( 'PDF' )
#
#
#histogram( values )
#moment_method( values )
#maximum_likelihood( values )
#plt.legend( loc = 'best' )
##plt.show()
