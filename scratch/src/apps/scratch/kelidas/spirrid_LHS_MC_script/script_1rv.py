from scipy.stats.distributions import norm, weibull_min # import normal distribution
import pylab as p  # import plotting tool
from numpy import vectorize, linspace, zeros_like, sign, sum as nsum, ones, \
                 corrcoef, sort, diff, array, exp, min, sqrt
from numpy.random import random
from time import clock

n_int = 10 # number of discretization points
n_k = 1 # number of random variables

# set the mean and standard deviation of the two random variables
xi_mean, xi_stdev = 0.019027, 0.0022891

# construct the normal distributions and get the methods
# for the evaluation of the probability density functions
s = 0.02
m = 10.
g_xi = weibull_min( m, scale = s ) #norm( loc = xi_mean, scale = xi_stdev )

# generate the grids for integration covering major part of the random domains
delta_xi = ( xi_mean + ( 4 * xi_stdev ) - xi_mean + ( 4 * xi_stdev ) ) / n_int
Theta_xi = linspace( xi_mean - ( 4 * xi_stdev ) + 0.5 * delta_xi, xi_mean + ( 4 * xi_stdev ) - 0.5 * delta_xi, n_int )
# LHS generate the grids for integration covering major part of the random domains
T_xi = g_xi.ppf( linspace( 0.5 / n_int, 1. - 0.5 / n_int, n_int ) ) #g_xi.ppf( linspace( 0.5 / n_int, 1. - 0.5 / n_int, n_int ) )
# MC generation
T_xi_MC = g_xi.rvs( n_int ** n_k )
#T_xi_MC = array( zip( *sorted( zip( random( n_int ** n_k ), g_xi.ppf( linspace( 0.5 / n_int, 1. - 0.5 / n_int, n_int ** n_k ) ) ) ) )[1] )

# grid spacing
d_xi = delta_xi#( Theta_xi[-1] - Theta_xi[0] ) / n_int

def Heaviside( x ):
    ''' Heaviside function '''
    return ( sign( x ) + 1.0 ) / 2.0

def q( e, xi ):
    ''' Response function of a single fiber '''
    return e * Heaviside( xi - e )

def q_ex( e ):
    ''' Response function of a single fiber '''
    return e * exp( -( e / 0.02 ) ** 10. )

# prepare the sequence of the control strains in a numpy array
e_arr = linspace( 0, .03, 1000 )
# define an array of the same size as e_arr
mu_q_arr = zeros_like( e_arr )


dG_xi = g_xi.pdf( Theta_xi ) * d_xi

dG_grid = dG_xi

def mu_q_e( e ):
    ''' Summation / integration  over the random domain '''
    q_e_grid = q( e, Theta_xi )
    q_dG_grid = q_e_grid * dG_grid # element by element product of two (m,m) arrays
    return nsum( q_dG_grid ) # nsum has been imported at line 3 from numpy 

mu_q_e_vct = vectorize( mu_q_e )

def sig_q_e( e ):
    ''' Summation / integration  over the random domain '''
    q_e_grid = q( e, Theta_xi )
    q_dG_grid = q_e_grid ** 2 * dG_grid
    return sqrt( nsum( q_dG_grid ) - mu_q_e( e ) ** 2 )

sig_q_e_vct = vectorize( sig_q_e )

start_time = clock()
mu_q_arr = mu_q_e_vct( e_arr )
print 'loop-less: elapsed time', clock() - start_time
sig_q_arr = sig_q_e_vct( e_arr )

def mu_q_e_LHS( e ):
    ''' Summation / integration  over the random domain '''
    q_e_grid = q( e, T_xi )
    return nsum( q_e_grid ) / n_int ** n_k

mu_q_e_LHS_vct = vectorize( mu_q_e_LHS )

def sig_q_e_LHS( e ):
    ''' Summation / integration  over the random domain '''
    q_e_grid = q( e, T_xi )
    q_dG_grid = q_e_grid ** 2 / n_int
    return sqrt( nsum( q_dG_grid ) - mu_q_e_LHS( e ) ** 2 )

sig_q_e_LHS_vct = vectorize( sig_q_e_LHS )

start_time = clock()
mu_q_arr_LHS = mu_q_e_LHS_vct( e_arr )
print 'loop-less: elapsed time', clock() - start_time
sig_q_arr_LHS = sig_q_e_LHS_vct( e_arr )

def mu_q_e_MC( e ):
    ''' Summation / integration  over the random domain '''
    q_e_grid = q( e, T_xi_MC )
    return nsum( q_e_grid ) / n_int ** n_k

mu_q_e_MC_vct = vectorize( mu_q_e_MC )

def sig_q_e_MC( e ):
    ''' Summation / integration  over the random domain '''
    q_e_grid = q( e, T_xi_MC )
    q_dG_grid = q_e_grid ** 2 / n_int
    return sqrt( nsum( q_dG_grid ) - mu_q_e_MC( e ) ** 2 )

sig_q_e_MC_vct = vectorize( sig_q_e_MC )

start_time = clock()
mu_q_arr_MC = mu_q_e_MC_vct( e_arr )
print 'loop-less: elapsed time', clock() - start_time
sig_q_arr_MC = sig_q_e_MC_vct( e_arr )


mu_q_ex = q_ex( e_arr )
sigm_ex = ( e_arr ** 2 * exp( -e_arr ** m * ( 1.0 / s ) ** m ) - e_arr ** 2 * exp( -e_arr ** m * ( 1.0 / s ) ** m ) ** 2 ) ** ( 1. / 2. )
#####################
# ERROR
#######################xx
print 'reg grid error', nsum( ( mu_q_arr - mu_q_ex ) ** 2 )
print 'LHS grid error', nsum( ( mu_q_arr_LHS - mu_q_ex ) ** 2 )
print 'MC error', nsum( ( mu_q_arr_MC - mu_q_ex ) ** 2 )


p.plot( e_arr, mu_q_ex, 'k-', linewidth = 2, label = 'exact' )
p.plot( e_arr, sigm_ex, 'k.', linewidth = 2, label = 'exact' )
p.plot( e_arr, mu_q_arr, 'b-', label = 'regular grid' )
p.plot( e_arr, sig_q_arr, 'b.-', label = 'regular grid' )
p.plot( e_arr, mu_q_arr_LHS, 'r-', label = 'LHS grid' )
p.plot( e_arr, sig_q_arr_LHS, 'r.-', label = 'LHS grid' )
p.plot( e_arr, mu_q_arr_MC, 'g-', label = 'MC' )
p.plot( e_arr, sig_q_arr_MC, 'g.-', label = 'MC' )
p.plot( Theta_xi, 0.003 * ones( n_int ), 'bo' )
p.plot( T_xi, 0.002 * ones( n_int ), 'ro' )
p.plot( T_xi_MC, 0.001 * ones( n_int ), 'go' )
p.legend( loc = 'best' )



#p.figure()
#p.plot( T_la[None, :] * ones( n_int )[:, None], T_xi[:, None] * ones( n_int )[None, :], 'ro', label = 'LHS grid' )
#p.plot( T_la_MC, T_xi_MC, 'gx', label = 'MC' )
#p.plot( Theta_la[None, :] * ones( n_int )[:, None], Theta_xi[:, None] * ones( n_int )[None, :], 'b.', label = 'regular grid' )


p.show()
