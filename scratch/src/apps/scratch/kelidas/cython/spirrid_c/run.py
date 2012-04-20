from scipy.stats.distributions import norm, weibull_min, uniform # import normal distribution
import pylab as p  # import plotting tool
from numpy import vectorize, linspace, zeros_like, sign, sum as nsum, ones, corrcoef, sort, diff, array, loadtxt
from numpy.random import random
from time import clock
from scipy.interpolate import interp1d

from scipy.weave import inline, converters


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
n_eps = 1000
e_arr = linspace( 0, 0.03, n_eps )


import spirrid

start_time = clock()
mu_q = spirrid.loop_mu_q_e( e_arr, d_la, d_xi, Theta_la, Theta_xi, g_la_pdf, g_xi_pdf )
t_cython = clock() - start_time
print 'loop-based: elapsed time', t_cython

print mu_q



C_code = '''
for( int i_eps = 0; i_eps < n_eps; i_eps++){
	double eps = *( e_arr + i_eps );
	double mu_q(0);
	double q(0);
	for( int i_lambd = 0; i_lambd < n_int; i_lambd++){
		double lambd = *( Theta_la + i_lambd );
	for( int i_xi = 0; i_xi < n_int; i_xi++){
		double xi = *( Theta_xi + i_xi );
		double pdf =  *( g_la_pdf + i_lambd)* *( g_xi_pdf + i_xi);

		    double eps_ = xi - ( eps ) /( 1 + lambd );
		    // Computation of the q( ... ) function
		    if ( eps_ < 0 || eps_ > xi ){
			q = 0.0;
		    }else{
			  q = ( eps ) /( 1 + lambd );
		    }
		// Store the values in the grid
		mu_q +=  q * pdf;
	};
	};
*(mu_q_arr + i_eps) = mu_q;
};
'''

mu_q_arr = zeros_like( e_arr )

arg_list = ['mu_q_arr', 'e_arr', 'Theta_la', 'Theta_xi', 'g_la_pdf', 'g_xi_pdf', 'n_eps', 'n_int']


c_params = {}
c_params['mu_q_arr'] = mu_q_arr
c_params['e_arr'] = e_arr
c_params['Theta_la'] = Theta_la
c_params['Theta_xi'] = Theta_xi
c_params['g_la_pdf'] = g_la_pdf
c_params['g_xi_pdf'] = g_xi_pdf
c_params['n_eps'] = n_eps
c_params['n_int'] = n_int

compiler = 'gcc'

compiler_verbose = 0

conv = converters.default


start_time = clock()
inline( C_code, arg_list, local_dict = c_params,
                                   type_converters = conv, compiler = compiler,
                                   verbose = compiler_verbose )
t_c = clock() - start_time
print 'loop-based: elapsed time', t_c


p.plot( e_arr, mu_q, label = 'Cython - %s s' % t_cython )
p.plot( e_arr, mu_q_arr * d_la * d_xi, label = 'C - %s s' % t_c )
p.legend()
p.show()
