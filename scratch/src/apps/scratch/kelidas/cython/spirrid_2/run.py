from scipy.stats.distributions import norm, weibull_min, uniform # import normal distribution
import pylab as p  # import plotting tool
from numpy import vectorize, linspace, zeros_like, sign, sum as nsum, ones, corrcoef, sort, diff, array, loadtxt
from numpy.random import random
from time import clock
from scipy.interpolate import interp1d


n_int = 10 # number of discretization points
n_k = 2 # number of random variables

# set the mean and standard deviation of the two random variables
la_mean, la_stdev = 0.0, 0.2
xi_mean, xi_stdev = 0.019027, 0.0022891

# construct the normal distributions and get the methods
# for the evaluation of the probability density functions
g_la = uniform( loc = la_mean, scale = la_stdev )
g_xi = weibull_min( 10., scale = 0.02 )

# generate the grids for integration covering major part of the random domains
Theta_la = linspace( la_mean + 0.5 * la_stdev / n_int, la_mean + la_stdev - 0.5 * la_stdev / n_int, n_int )
delta_xi = ( xi_mean + ( 4 * xi_stdev ) - xi_mean + ( 4 * xi_stdev ) ) / n_int
Theta_xi = linspace( xi_mean - ( 4 * xi_stdev ) + 0.5 * delta_xi, xi_mean + ( 4 * xi_stdev ) - 0.5 * delta_xi, n_int )
# LHS generate the grids for integration covering major part of the random domains
T_la = g_la.ppf( linspace( 0.5 / n_int, 1. - 0.5 / n_int, n_int ) )
T_xi = g_xi.ppf( linspace( 0.5 / n_int, 1. - 0.5 / n_int, n_int ) )
# MC generation
T_la_MC = g_la.rvs( n_int ** n_k )
T_xi_MC = g_xi.rvs( n_int ** n_k )
#T_la_MC = array( zip( *sorted( zip( random( n_int ** n_k ), g_la.ppf( linspace( 0.5 / n_int, 1. - 0.5 / n_int, n_int ** n_k ) ) ) ) )[1] )
#T_xi_MC = array( zip( *sorted( zip( random( n_int ** n_k ), g_xi.ppf( linspace( 0.5 / n_int, 1. - 0.5 / n_int, n_int ** n_k ) ) ) ) )[1] )
print diff( sort( g_la.cdf( g_la.rvs( 10 ) ) ) )
print 'MC - correlation between la and xi', corrcoef( T_la_MC, T_xi_MC )

# grid spacing
d_la = la_stdev / n_int#( Theta_la[-1] - Theta_la[0] ) / n_int
d_xi = delta_xi#( Theta_xi[-1] - Theta_xi[0] ) / n_int 

g_la_pdf = g_la.pdf( Theta_la )
g_xi_pdf = g_xi.pdf( Theta_xi )

# prepare the sequence of the control strains in a numpy array
e_arr = linspace( 0, 0.03, 1000 )


import spirrid

start_time = clock()
mu_q = spirrid.loop_mu_q_e( e_arr, d_la, d_xi, Theta_la, Theta_xi, g_la_pdf, g_xi_pdf )
print 'loop-based: elapsed time', clock() - start_time
print mu_q


p.plot( e_arr, mu_q )
p.show()
