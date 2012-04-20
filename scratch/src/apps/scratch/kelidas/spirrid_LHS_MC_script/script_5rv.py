from scipy.stats.distributions import norm, uniform # import normal distribution
import pylab as p  # import plotting tool
from numpy import vectorize, linspace, zeros_like, sign, sum as nsum, ones, \
                 corrcoef, sort, diff, array, pi, trapz, sqrt, loadtxt
from numpy.random import random
from time import clock
from scipy.interpolate import interp1d

n_int = 4 # number of discretization points
n_k = 5 # number of random variables

D = 26 * 1.0e-6 # m
A = ( D / 2.0 ) ** 2 * pi

# set the mean and standard deviation of the two random variables
la_mean, la_stdev = 0.0, 0.2
xi_mean, xi_stdev = 0.019027, 0.0022891
E_mean, E_stdev = 70.0e+9, 15.0e+9
th_mean, th_stdev = 0.0, 0.01
A_mean, A_stdev = A * 0.3, 0.7 * A
print A_mean, A_mean + A_stdev

# construct the normal distributions and get the methods
# for the evaluation of the probability density functions
g_la = uniform( loc = la_mean, scale = la_stdev )
g_xi = norm( loc = xi_mean, scale = xi_stdev )
g_E = uniform( loc = E_mean, scale = E_stdev )
g_th = uniform( loc = th_mean, scale = th_stdev )
g_A = uniform( loc = A_mean, scale = A_stdev )

# generate the grids for integration covering major part of the random domains
Theta_la = linspace( la_mean + 0.5 * la_stdev / n_int, la_mean + la_stdev - 0.5 * la_stdev / n_int, n_int )
delta = ( xi_mean + ( 4 * xi_stdev ) - xi_mean + ( 4 * xi_stdev ) ) / n_int
Theta_xi = linspace( xi_mean - ( 4 * xi_stdev ) + 0.5 * delta, xi_mean + ( 4 * xi_stdev ) - 0.5 * delta, n_int )
Theta_E = linspace( E_mean + 0.5 * E_stdev / n_int, E_mean + E_stdev - 0.5 * E_stdev / n_int, n_int )
Theta_th = linspace( th_mean + 0.5 * th_stdev / n_int, th_mean + th_stdev - 0.5 * th_stdev / n_int, n_int )
Theta_A = linspace( A_mean + 0.5 * A_stdev / n_int, A_mean + A_stdev - 0.5 * A_stdev / n_int, n_int )
# LHS generate the grids for integration covering major part of the random domains
T_la = g_la.ppf( linspace( 0.5 / n_int, 1. - 0.5 / n_int, n_int ) )
T_xi = g_xi.ppf( linspace( 0.5 / n_int, 1. - 0.5 / n_int, n_int ) )
T_E = g_E.ppf( linspace( 0.5 / n_int, 1. - 0.5 / n_int, n_int ) )
T_th = g_th.ppf( linspace( 0.5 / n_int, 1. - 0.5 / n_int, n_int ) )
T_A = g_A.ppf( linspace( 0.5 / n_int, 1. - 0.5 / n_int, n_int ) )
# MC generation
T_la_MC = g_la.rvs( n_int ** n_k )
T_xi_MC = g_xi.rvs( n_int ** n_k )
T_E_MC = g_E.rvs( n_int ** n_k )
T_th_MC = g_th.rvs( n_int ** n_k )
T_A_MC = g_A.rvs( n_int ** n_k )
#T_la_MC = array( zip( *sorted( zip( random( n_int ** n_k ), g_la.ppf( linspace( 0.5 / n_int, 1. - 0.5 / n_int, n_int ** n_k ) ) ) ) )[1] )
#T_xi_MC = array( zip( *sorted( zip( random( n_int ** n_k ), g_xi.ppf( linspace( 0.5 / n_int, 1. - 0.5 / n_int, n_int ** n_k ) ) ) ) )[1] )
#T_E_MC = array( zip( *sorted( zip( random( n_int ** n_k ), g_E.ppf( linspace( 0.5 / n_int, 1. - 0.5 / n_int, n_int ** n_k ) ) ) ) )[1] )
#T_th_MC = array( zip( *sorted( zip( random( n_int ** n_k ), g_th.ppf( linspace( 0.5 / n_int, 1. - 0.5 / n_int, n_int ** n_k ) ) ) ) )[1] )
#T_A_MC = array( zip( *sorted( zip( random( n_int ** n_k ), g_A.ppf( linspace( 0.5 / n_int, 1. - 0.5 / n_int, n_int ** n_k ) ) ) ) )[1] )
print 'MC - correlation between la and xi', corrcoef( T_la_MC, T_xi_MC )

# grid spacing
d_la = la_stdev / n_int #( Theta_la[-1] - Theta_la[0] ) / n_int
d_xi = delta #( Theta_xi[-1] - Theta_xi[0] ) / n_int
d_E = E_stdev / n_int #( Theta_E[-1] - Theta_E[0] ) / n_int
d_th = th_stdev / n_int #( Theta_th[-1] - Theta_th[0] ) / n_int
d_A = A_stdev / n_int #( Theta_A[-1] - Theta_A[0] ) / n_int
print d_E, d_th, d_A, d_xi, d_la

def Heaviside( x ):
    ''' Heaviside function '''
    return ( sign( x ) + 1.0 ) / 2.0

def q( eps, lambd, xi, E_mod, theta, A ):
    eps_ = ( eps - theta * ( 1 + lambd ) ) / ( ( 1 + theta ) * ( 1 + lambd ) )
    eps_ *= Heaviside( eps_ )
    eps_grid = eps_ * Heaviside( xi - eps_ )
    q_grid = E_mod * A * eps_grid
    return q_grid

def q_ex( e ):
    data = loadtxt( 'mu_q_5rv_30.txt', delimiter = ',' )
    x, y = data[:, 0], data[:, 1]
    f = interp1d( x, y, kind = 'linear' )
    return f( e )

def sig_ex( e ):
    data = loadtxt( 'sig_5rv_30.txt', delimiter = ',' )
    x, y = data[:, 0], data[:, 1]
    f = interp1d( x, y, kind = 'linear' )
    return f( e )

# prepare the sequence of the control strains in a numpy array
e_arr = linspace( 0, 0.05, 1000 )
# define an array of the same size as e_arr
mu_q_arr = zeros_like( e_arr )

def loop_mu_q_e():
    # loop over the control variable (strain)
    for i, e in enumerate( e_arr ):
        mu_q_e = 0.0  # interim variable for the summed stress contributions
        for la in Theta_la: # loop over lambda range (array of values)
            for xi in Theta_xi: # loop over xi range (array of values)
                for E in Theta_E:
                    for th in Theta_th:
                        for A in Theta_A:
                            dG = ( g_la.pdf( la ) * g_xi.pdf( xi ) *
                                   g_E.pdf( E ) * g_th.pdf( th ) *
                                    g_A.pdf( A ) *
                                   d_la * d_xi * d_E * d_th * d_A )
                            mu_q_e += q( e, la, xi, E, th, A ) * dG
        mu_q_arr[ i ] = mu_q_e
    return mu_q_arr

start_time = clock()
#mu_q_arr_loop = loop_mu_q_e()

print 'loop-based: elapsed time', clock() - start_time

dG_la = g_la.pdf( Theta_la ) * d_la
dG_xi = g_xi.pdf( Theta_xi ) * d_xi
dG_E = g_E.pdf( Theta_E ) * d_E
dG_th = g_th.pdf( Theta_th ) * d_th
dG_A = g_A.pdf( Theta_A ) * d_A

print nsum( dG_la ), nsum( dG_xi ), nsum( dG_E ), nsum( dG_th ), nsum( dG_A )

dG_grid = ( dG_la[:, None, None, None, None] * dG_xi[None, :, None, None, None] *
           dG_E[None, None, :, None, None] * dG_th[None, None, None, :, None] *
            dG_A[None, None, None, None, :] )


def mu_q_e( e ):
    ''' Summation / integration  over the random domain '''
    q_e_grid = q( e, Theta_la[:, None, None, None, None], Theta_xi[None, :, None, None, None],
                   Theta_E[None, None, :, None, None], Theta_th[None, None, None, :, None],
                   Theta_A[None, None, None, None, :] )
    q_dG_grid = q_e_grid * dG_grid # element by element product of two (m,m) arrays
    return nsum( q_dG_grid ) # nsum has been imported at line 3 from numpy 

mu_q_e_vct = vectorize( mu_q_e )

def sig_q_e( e ):
    ''' Summation / integration  over the random domain '''
    q_e_grid = q( e, Theta_la[:, None, None, None, None], Theta_xi[None, :, None, None, None],
                   Theta_E[None, None, :, None, None], Theta_th[None, None, None, :, None],
                   Theta_A[None, None, None, None, :] )
    q_dG_grid = q_e_grid ** 2 * dG_grid
    return sqrt( nsum( q_dG_grid ) - mu_q_e( e ) ** 2 )

sig_q_e_vct = vectorize( sig_q_e )


start_time = clock()
mu_q_arr = mu_q_e_vct( e_arr )
print 'loop-less: elapsed time', clock() - start_time
sig_q_arr = sig_q_e_vct( e_arr )

def mu_q_e_LHS( e ):
    ''' Summation / integration  over the random domain '''
    q_e_grid = q( e, T_la[:, None, None, None, None], T_xi[None, :, None, None, None],
                   T_E[None, None, :, None, None], T_th[None, None, None, :, None],
                   T_A[None, None, None, None, :] )
    return nsum( q_e_grid ) / n_int ** n_k

mu_q_e_LHS_vct = vectorize( mu_q_e_LHS )

def sig_q_e_LHS( e ):
    ''' Summation / integration  over the random domain '''
    q_e_grid = q( e, T_la[:, None, None, None, None], T_xi[None, :, None, None, None],
                   T_E[None, None, :, None, None], T_th[None, None, None, :, None],
                   T_A[None, None, None, None, :] )
    q_dG_grid = q_e_grid ** 2 / n_int ** n_k
    return sqrt( nsum( q_dG_grid ) - mu_q_e_LHS( e ) ** 2 )

sig_q_e_LHS_vct = vectorize( sig_q_e_LHS )


start_time = clock()
mu_q_arr_LHS = mu_q_e_LHS_vct( e_arr )
print 'loop-less: elapsed time', clock() - start_time
sig_q_arr_LHS = sig_q_e_LHS_vct( e_arr )


def mu_q_e_MC( e ):
    ''' Summation / integration  over the random domain '''
    q_e_grid = q( e, T_la_MC, T_xi_MC, T_E_MC, T_th_MC, T_A_MC )
    return nsum( q_e_grid ) / n_int ** n_k

mu_q_e_MC_vct = vectorize( mu_q_e_MC )

def sig_q_e_MC( e ):
    ''' Summation / integration  over the random domain '''
    q_e_grid = q( e, T_la_MC, T_xi_MC, T_E_MC, T_th_MC, T_A_MC )
    q_dG_grid = q_e_grid ** 2 / n_int ** n_k
    return sqrt( nsum( q_dG_grid ) - mu_q_e_MC( e ) ** 2 )

sig_q_e_MC_vct = vectorize( sig_q_e_MC )


start_time = clock()
mu_q_arr_MC = mu_q_e_MC_vct( e_arr )
print 'loop-less: elapsed time', clock() - start_time
sig_q_arr_MC = sig_q_e_MC_vct( e_arr )




mu_q_ex = q_ex( e_arr )
sig_q_ex = sig_ex( e_arr )




p.plot( e_arr, mu_q_ex, 'k-', linewidth = 2, label = 'n_int = 30' )
p.plot( e_arr, sig_q_ex, 'k-', linewidth = 2, label = 'n_int = 30' )
#p.plot( e_arr, mu_q_arr_loop, 'c-x', label = 'regular grid loop' )
p.plot( e_arr, mu_q_arr, 'b-', label = 'regular grid' )
p.plot( e_arr, sig_q_arr, 'b.-', label = 'regular grid' )
p.plot( e_arr, mu_q_arr_LHS, 'r-', label = 'LHS grid' )
p.plot( e_arr, sig_q_arr_LHS, 'r.-', label = 'LHS grid' )
p.plot( e_arr, mu_q_arr_MC, 'g-', label = 'MC' )
p.plot( e_arr, sig_q_arr_MC, 'g.-', label = 'MC' )
p.plot( Theta_xi, 0.05 * ones( n_int ), 'bo' )
p.plot( 0, 0, 'b|', label = '$\\xi(1+\lambda)(1+\\theta)+\\theta(1+\lambda)$' )
p.plot( ( Theta_xi[None, None, :] *
           ( 1 + Theta_th[None, :, None] ) *
            ( 1 + Theta_la[:, None, None] ) +
         Theta_th[None, :, None] *
         ( 1 + Theta_la[:, None, None] ) ).flatten(),
        0.04 * ones( n_int ** 3 ), 'b|' )
p.plot( T_xi, 0.1 * ones( n_int ), 'ro' )
p.plot( ( T_xi[None, None, :] *
           ( 1 + T_th[None, :, None] ) *
            ( 1 + T_la[:, None, None] ) +
         T_th[None, :, None] *
         ( 1 + T_la[:, None, None] ) ).flatten(),
        0.09 * ones( n_int ** 3 ), 'r|' )
p.plot( T_xi_MC, 0.15 * ones( n_int ** n_k ), 'go' )
p.plot( ( T_xi_MC *
           ( 1 + T_th_MC ) *
            ( 1 + T_la_MC ) +
         T_th_MC *
         ( 1 + T_la_MC ) ),
        0.14 * ones( n_int ** 5 ), 'g|' )
p.legend()

p.figure()
p.plot( T_la[None, :] * ones( n_int )[:, None], T_xi[:, None] * ones( n_int )[None, :], 'ro', label = 'LHS grid' )
p.plot( T_la_MC, T_xi_MC, 'gx', label = 'MC' )
p.plot( Theta_la[None, :] * ones( n_int )[:, None], Theta_xi[:, None] * ones( n_int )[None, :], 'b.', label = 'regular grid' )


#p.figure()
#p.plot( Theta_xi, ones( n_int ) * 1, 'b,' )
#p.plot( T_xi, ones( n_int ) * 2, 'r,' )
#p.plot( T_xi_MC, ones( n_int ** n_k ) * 3, 'g,' )

p.show()
