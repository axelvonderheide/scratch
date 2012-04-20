from scipy.stats.distributions import norm # import normal distribution
import pylab as p  # import plotting tool
from numpy import vectorize, linspace, zeros_like, sign, sum as nsum

n_int = 20 # number of sampling points

# set the mean and the range of the sampling
c_mean, c_stdev = 1.0, 0.2
x_mean, x_stdev = 10.0, 2.0

# generate arrays of sampling values
c_arr = linspace( c_mean - ( 4 * c_stdev ), c_mean + ( 4 * c_stdev ), n_int )
x_arr = linspace( x_mean - ( 4 * x_stdev ), x_mean + ( 4 * x_stdev ), n_int )

# grid distances
dc = ( c_arr[-1] - c_arr[0] ) / n_int
dx = ( x_arr[-1] - x_arr[0] ) / n_int

# construct the normal distribution and get the method
# for the evaluation of the cumulative probability
pdf_c = norm( loc = c_mean, scale = c_stdev ).pdf
pdf_x = norm( loc = x_mean, scale = x_stdev ).pdf

def Heaviside( x ):
    ''' Heaviside function '''
    #@TODO: same as definition
    return ( sign( x ) + 1.0 ) / 2.

def q( e, c, x ):
    ''' Response function of a single fiber '''
    return c * e * Heaviside( x - e )

# prepare the sequence of the control strains
# evaluate the response for an array of values of the control variable
e_arr = linspace( 0, 20, 100 )
# define an array of the same size as e_arr
mu_q_arr = zeros_like( e_arr )

def loop_mu_q_e():
    # loop over the control variable (strain)
    for i, e in enumerate( e_arr ):
        mu_q_e = 0.0  # interim variable for the summed stress contributions
        for c in c_arr: # loop over lambda range (array of values)
            for x in x_arr: # loop over xi range (array of values)
                dG = pdf_c( c ) * pdf_x( x ) * dc * dx
                mu_q_e += q( e, c, x ) * dG
        mu_q_arr[ i ] = mu_q_e

dG_c = pdf_c( c_arr ) * dc
dG_x = pdf_x( x_arr ) * dx

dG_grid = dG_c[:, None] * dG_x[None, :]

def mu_q_e( e ):
    ''' Summation / integration  over the random domain '''
    q_e_grid = q( e, c_arr[:, None], x_arr[None, :] )
    q_dG_grid = q_e_grid * dG_grid
    return nsum( q_dG_grid )

mu_q_e_vct = vectorize( mu_q_e )

mu_q_arr = mu_q_e_vct( e_arr )
loop_mu_q_e()
p.plot( e_arr, mu_q_arr )
p.show()
