from scipy.stats.distributions import norm, weibull_min, uniform # import normal distribution
import pylab as p  # import plotting tool
from numpy import vectorize, linspace, zeros_like, sign, sum as nsum, ones, \
                 corrcoef, sort, diff, array, loadtxt, sqrt
from numpy.random import random
from time import clock
from scipy.interpolate import interp1d


n_int = 5 # number of discretization points
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

def Heaviside( x ):
    ''' Heaviside function '''
    return ( sign( x ) + 1.0 ) / 2.0

def q( e, la, xi ):
    ''' Response function of a single fiber '''
    return  e / ( 1 + la ) * Heaviside( xi - e / ( 1 + la ) )

def q_ex( e ):
    data = loadtxt( '2rv_maple.txt', delimiter = ',' )
    x, y = data[:, 0], data[:, 1]
    f = interp1d( x, y, kind = 'linear' )
    return f( e )

def sig_ex( e ):
    data = loadtxt( '2rv_maple_sig.txt', delimiter = ',' )
    x, y = data[:, 0], data[:, 1]
    f = interp1d( x, y, kind = 'linear' )
    return f( e )

# prepare the sequence of the control strains in a numpy array
e_arr = linspace( 0, 0.03, 1000 )
# define an array of the same size as e_arr
mu_q_arr = zeros_like( e_arr )

def loop_mu_q_e():
    # loop over the control variable (strain)
    for i, e in enumerate( e_arr ):
        mu_q_e = 0.0  # interim variable for the summed stress contributions
        for la in Theta_la: # loop over lambda range (array of values)
            for xi in Theta_xi: # loop over xi range (array of values)
                dG = g_la.pdf( la ) * g_xi.pdf( xi ) * d_la * d_xi
                mu_q_e += q( e, la, xi ) * dG
        mu_q_arr[ i ] = mu_q_e
    return mu_q_arr

start_time = clock()
#mu_q_loop = loop_mu_q_e()

print 'loop-based: elapsed time', clock() - start_time

dG_la = g_la.pdf( Theta_la ) * d_la
dG_xi = g_xi.pdf( Theta_xi ) * d_xi

dG_grid = dG_la[:, None] * dG_xi[None, :]

def mu_q_e( e ):
    ''' Summation / integration  over the random domain '''
    q_e_grid = q( e, Theta_la[:, None], Theta_xi[None, :] )
    q_dG_grid = q_e_grid * dG_grid # element by element product of two (m,m) arrays
    return nsum( q_dG_grid ) # nsum has been imported at line 3 from numpy 

mu_q_e_vct = vectorize( mu_q_e )

def sig_q_e( e ):
    ''' Summation / integration  over the random domain '''
    q_e_grid = q( e, Theta_la[:, None], Theta_xi[None, :] )
    q_dG_grid = q_e_grid ** 2 * dG_grid
    return sqrt( nsum( q_dG_grid ) - mu_q_e( e ) ** 2 )

sig_q_e_vct = vectorize( sig_q_e )

start_time = clock()
mu_q_arr = mu_q_e_vct( e_arr )
print 'loop-less: elapsed time', clock() - start_time
sig_q_arr = sig_q_e_vct( e_arr )


def mu_q_e_LHS( e ):
    ''' Summation / integration  over the random domain '''
    q_e_grid = q( e, T_la[:, None], T_xi[None, :] )
    return nsum( q_e_grid ) / n_int ** n_k

mu_q_e_LHS_vct = vectorize( mu_q_e_LHS )

def sig_q_e_LHS( e ):
    ''' Summation / integration  over the random domain '''
    q_e_grid = q( e, T_la[:, None], T_xi[None, :] )
    q_dG_grid = q_e_grid ** 2 / n_int ** n_k
    return sqrt( nsum( q_dG_grid ) - mu_q_e_LHS( e ) ** 2 )

sig_q_e_LHS_vct = vectorize( sig_q_e_LHS )

start_time = clock()
mu_q_arr_LHS = mu_q_e_LHS_vct( e_arr )
print 'loop-less: elapsed time', clock() - start_time
sig_q_arr_LHS = sig_q_e_LHS_vct( e_arr )


def mu_q_e_MC( e ):
    ''' Summation / integration  over the random domain '''
    q_e_grid = q( e, T_la_MC, T_xi_MC )
    return nsum( q_e_grid ) / n_int ** n_k

mu_q_e_MC_vct = vectorize( mu_q_e_MC )

def sig_q_e_MC( e ):
    ''' Summation / integration  over the random domain '''
    q_e_grid = q( e, T_la_MC, T_xi_MC )
    q_dG_grid = q_e_grid ** 2 / n_int ** n_k
    return sqrt( nsum( q_dG_grid ) - mu_q_e_MC( e ) ** 2 )

sig_q_e_MC_vct = vectorize( sig_q_e_MC )

start_time = clock()
mu_q_arr_MC = mu_q_e_MC_vct( e_arr )
print 'loop-less: elapsed time', clock() - start_time
sig_q_arr_MC = sig_q_e_MC_vct( e_arr )


mu_q_ex = q_ex( e_arr )
sig_q_ex = sig_ex( e_arr )

#####################
# ERROR
#######################xx
print 'reg grid error', nsum( ( mu_q_arr - mu_q_ex ) ** 2 )
print 'LHS grid error', nsum( ( mu_q_arr_LHS - mu_q_ex ) ** 2 )
print 'MC error', nsum( ( mu_q_arr_MC - mu_q_ex ) ** 2 )

#p.plot( e_arr, mu_q_loop, 'r.' )
p.plot( e_arr, mu_q_ex, 'k-', label = 'maple' )
p.plot( e_arr, sig_q_ex, 'k-', label = 'maple' )
p.plot( e_arr, mu_q_arr, 'b-', label = 'regular grid' )
p.plot( e_arr, sig_q_arr, 'b.-', label = 'regular grid' )
p.plot( e_arr, mu_q_arr_LHS, 'r-', label = 'LHS grid' )
p.plot( e_arr, sig_q_arr_LHS, 'r.-', label = 'LHS grid' )
p.plot( e_arr, mu_q_arr_MC, 'g-', label = 'MC' )
p.plot( e_arr, sig_q_arr_MC, 'g.-', label = 'MC' )
p.plot( Theta_xi, 0.003 * ones( n_int ), 'bo' )
p.plot( Theta_xi * ( 1 + Theta_la ), .0028 * ones( n_int ), 'b|', label = '$\\xi(1+\lambda)$' )
p.plot( Theta_xi[None, :] * ( 1 + Theta_la[:, None] ), 0.0028 * ones( n_int ), 'b|' )
p.plot( T_xi, 0.002 * ones( n_int ), 'ro' )
p.plot( T_xi * ( 1 + T_la ), 0.0018 * ones( n_int ), 'r|', label = '$\\xi(1+\lambda)$' )
p.plot( T_xi[None, :] * ( 1 + T_la[:, None] ), 0.0018 * ones( n_int ), 'r|' )
p.plot( T_xi_MC, 0.001 * ones( n_int ** 2 ), 'go' )
p.plot( 0, 0, 'g|', label = '$\\xi(1+\lambda)$' )
p.plot( T_xi_MC * ( 1 + T_la_MC ), 0.0008 * ones( n_int ** 2 ), 'g|' )
p.legend()

p.figure()
p.plot( T_la[None, :] * ones( n_int )[:, None], T_xi[:, None] * ones( n_int )[None, :], 'ro', label = 'LHS grid' )
p.plot( T_la_MC, T_xi_MC, 'gx', label = 'MC' )
p.plot( Theta_la[None, :] * ones( n_int )[:, None], Theta_xi[:, None] * ones( n_int )[None, :], 'b.', label = 'regular grid' )


p.show()
