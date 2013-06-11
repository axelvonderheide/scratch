
from scipy.stats.distributions import norm
from numpy import linspace, sign, sum as nsum, vectorize

# construct the dG multidimensional arrayp
la_mean, la_range = 1.0, 0.2
xi_mean, xi_range = 10.0, 4.0

cdf_la = norm( loc = la_mean, scale = 0.2 * la_mean ).cdf
cdf_xi = norm( loc = xi_mean, scale = 0.2 * xi_mean ).cdf

n_int = 30.0
range_la = linspace( la_mean - la_range, la_mean + la_range, n_int )
range_xi = linspace( xi_mean - xi_range, xi_mean + xi_range, n_int )

dG_la = cdf_la( range_la ) * la_range * 2.0 / n_int
dG_xi = cdf_xi( range_xi ) * xi_range * 2.0 / n_int

dG_grid = dG_la[:, None] * dG_xi[None, :]

def H( x ):
    ''' Heaviside function '''
    return ( sign( x ) + 1.0 ) / 2.

def q( e, la, xi ):
    ''' Response function of a single fiber '''
    return e / ( 1 + la ) * H( xi - e / ( 1 + la ) )

def mu_q_e( e ):
    ''' Summation / integration  over the random domain '''
    q_e_grid = q( e, range_la[:, None], range_xi[None, :] )
    q_dG_grid = q_e_grid * dG_grid
    return nsum( q_dG_grid )

mu_q_fn = vectorize( mu_q_e )

# evaluate the response for an array of values of the control variable
e_arr = linspace( 0, 40, 100 )
mu_q_arr = mu_q_fn( e_arr )

####################################3

def plot_curves( p, e_arr ):
    for la in range_la:
        for xi in range_xi:
            p.plot( e_arr, q( e_arr, la, xi ), color = 'gray' )

import pylab as p
plot_curves( p, e_arr )
p.plot( e_arr, mu_q_fn( e_arr ), linewidth = 3 )
p.show()
