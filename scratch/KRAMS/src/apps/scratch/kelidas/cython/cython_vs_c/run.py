
from c_code import spirrid_c
from numpy import linspace, zeros_like
from scipy.stats.distributions import weibull_min, uniform # import normal distribution
from time import clock
import pylab as p  # import plotting tool

import pyximport; pyximport.install()
import spirrid_cython as spirrid

n_int = 1000 # number of discretization points
n_k = 2 # number of random variables 

# set the mean and standard deviation of the two random variables
la_mean, la_stdev = 0.0, 0.2
xi_mean, xi_stdev = 0.019027, 0.0022891

# construct the normal distributions and get the methods
# for the evaluation of the probability density functions
g_la = uniform(loc = la_mean, scale = la_stdev)
g_xi = weibull_min(10., scale = 0.02)

# generate the grids for integration covering major part of the random domains
theta_arr = linspace(-(1.0 - 1.0 / n_int), 1.0 - 1.0 / n_int, n_int)
Theta_la = la_mean + 4 * la_stdev * theta_arr
Theta_xi = xi_mean + 4 * xi_stdev * theta_arr

# get the size of the integration cell
d_la = (8 * la_stdev) / n_int
d_xi = (8 * xi_stdev) / n_int

g_la_pdf = g_la.pdf(Theta_la)
g_xi_pdf = g_xi.pdf(Theta_xi)

# prepare the sequence of the control strains in a numpy array
n_eps = 1000
e_arr = linspace(0, 0.03, n_eps)

mu_q_arr = zeros_like(e_arr)

arg_list = ['mu_q_arr', 'e_arr', 'Theta_la', 'Theta_xi', 'g_la_pdf', 'g_xi_pdf', 'n_eps', 'n_int', 'd_la', 'd_xi']


c_params = {}
c_params['mu_q_arr'] = mu_q_arr
c_params['e_arr'] = e_arr
c_params['Theta_la'] = Theta_la
c_params['Theta_xi'] = Theta_xi
c_params['g_la_pdf'] = g_la_pdf
c_params['g_xi_pdf'] = g_xi_pdf
c_params['n_eps'] = n_eps
c_params['n_int'] = n_int
c_params['d_la'] = d_la
c_params['d_xi'] = d_xi

#===============================================================================
# RUN C-code and Cython
#===============================================================================

# run C-code
start_time = clock()
spirrid_c(arg_list, c_params)
t_c = clock() - start_time
print 'c_code', t_c

# run Cython
start_time = clock()
mu_q = spirrid.loop_mu_q_e(e_arr, d_la, d_xi, Theta_la, Theta_xi, g_la_pdf, g_xi_pdf)
t_cython = clock() - start_time
print 'cython', t_cython

# plot results
p.plot(e_arr, mu_q_arr, 'r,', label = 'C - %s s' % t_c)
p.plot(e_arr, mu_q, 'b-', label = 'Cython - %s s' % t_cython)
p.legend()
p.show()
