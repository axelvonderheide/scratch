'''
@author: axel
'''
from math import pi as Pi
import numpy as np
from scipy.stats import norm
from stats.pdistrib.sin2x_distr import sin2x
from quaducom.resp_func.cb_short_fiber import CBShortFiber
from scipy.stats import exponweib
from matplotlib import pyplot as plt

#set parameters
Le = 4.5
cbsf = CBShortFiber()
max = 0.005
nxz = 1000
U = np.linspace( 0, max, nxz )
anzahl_der_realisationen = 25
nx = 50

length = 100
width = 30
height = 30
fiber_length = 9
r = 0.075
f = 0.86
def H( x ):
    return np.sign( np.sign( x ) + 1. )

start = 0.0001
end = 0.99999
dp = 3000
############################################


x_array = np.linspace( start, end, dp )

def pdf_tau( x ) :
    return  exponweib( 7.38493171518, 2.01255297859 ).pdf( x )


def ppf_tau ( x ):
    return exponweib( 7.38493171518, 2.01255297859 ).ppf( x )

def pdf_phi ( x ):
    return sin2x.pdf( x )

def ppf_phi ( x ):
    return sin2x.ppf( x )


def pdf_le ( x ):
    return 1. / Le * ( x ** 0 )


def ppf_le ( x ):
    return x * Le




#Realisationen der Zufallsparameter

real_tau_or = np.linspace( ppf_tau( start ), ppf_tau( end ), nx )
real_tau = real_tau_or.reshape( 1, 1, len( real_tau_or ) )

real_phi_or = np.linspace( ppf_phi( start ), ppf_phi( end ), nx )
real_phi = real_phi_or.reshape( 1, len( real_phi_or ), 1 )

#real_f_or = np.linspace( ppf_f( start ), ppf_f( end ), nx )
#real_f = real_f_or.reshape( 1, len( real_f_or ), 1, 1 )

real_le_or = np.linspace( ppf_le( start ), ppf_le( end ), nx )
real_le = real_le_or.reshape( len( real_le_or ), 1, 1 )


x_tau = np.linspace( ppf_tau( start ), ppf_tau( end ), nx )
x_phi = np.linspace( ppf_phi( start ), ppf_phi( end ), nx )
#x_f = np.linspace( ppf_f( start ), ppf_f( end ), nx )
x_le = np.linspace( ppf_le( start ), ppf_le( end ), nx )
#Distribution

z_tau = pdf_tau( real_tau )
z_phi = pdf_phi( real_phi )
#z_f = pdf_f( real_f )
z_le = pdf_le( real_le )
test_x = np.linspace( 0, 10, 30 )


pdf_matrix = z_tau * z_phi * z_le


# Integration
mean_matrix_list = []
variance_list = []
for w in U:
    pullouts = cbsf( w, real_tau, 0, 0.15, 200000, real_le, real_phi, f, 0, 0, 1e15 )
    mean_to_be_integrated = pullouts * pdf_matrix
    #print to_be_integrated
    step_one = np.trapz( mean_to_be_integrated, x_tau )
    step_two = np.trapz( step_one, x_phi )
    step_four = np.trapz( step_two, x_le )
    mean_matrix_list.append( step_four )

    stdev_to_be_integrated = ( pullouts - step_four ) ** 2 * pdf_matrix
    step_one1 = np.trapz( stdev_to_be_integrated, x_tau )
    step_two2 = np.trapz( step_one1, x_phi )
    step_four4 = np.trapz( step_two2, x_le )
    variance_list.append( step_four4 )


mean_array = np.array( mean_matrix_list )
var_fiber_array = np.array( variance_list )
stdev_array = var_fiber_array ** 0.5
lower_bound_scatter_band_pre = mean_array - stdev_array
lower_bound_scatter_band = H( lower_bound_scatter_band_pre ) * lower_bound_scatter_band_pre
upper_bound_scatter_band = mean_array + stdev_array
lower_twice_stdev = mean_array - stdev_array * 2
upper_twice_stdev = mean_array + stdev_array * 2
lower_triple_stdev = mean_array - stdev_array * 3
upper_triple_stdev = mean_array + stdev_array * 3


#plt.plot( U, lower_twice_stdev, 'k--', U, upper_twice_stdev, 'k--' )
#plt.plot( U, lower_triple_stdev, 'k--', U, upper_triple_stdev, 'k--' )
plt.title( 'Kraft-Rissoeffnungskurve', fontsize = 16 )
plt.xlabel( 'Verschiebung $w$ in [mm]', fontsize = 16 )
plt.ylabel( 'Kraft in [N]', fontsize = 16 )
#one_stdev = plt.fill_between( U, lower_triple_stdev , upper_triple_stdev  , color = '0.9', label = '99 Umgebung' )
#twice_stdev = plt.fill_between( U, lower_twice_stdev, upper_twice_stdev , color = '0.6' )
#triple_stdev = plt.fill_between( U, lower_bound_scatter_band, upper_bound_scatter_band , color = '0.3' )
#plt.legend( ( one_stdev, twice_stdev, triple_stdev ), ( '99 Umgebung', '95 Umgebung', '68 Umgebung' ), 'lower right' )
w = np.linspace( 0, max, nxz )
U = U.reshape( nxz, 1 )
for i in range( anzahl_der_realisationen ):

    tau_realisationen = ppf_tau( np.random.rand( 1 ) )
    phi_realisationen = ppf_phi( np.random.rand( 1 ) )
    #f_realisationen = ppf_phi( np.random.rand( 1 ) )
    phi_realisationen = ppf_phi( np.random.rand( 1 ) )
    le_realisationen = ppf_le( np.random.rand( 1 ) )
    pullouts = cbsf( U, tau_realisationen, 0, 0.15, 200000, le_realisationen, phi_realisationen, f, 0, 0, 1e15 )
    results = np.sum( pullouts , axis = 1 )
    plt.plot( w, results , 'b' )
z = plt.plot( w, results , 'b' )
a = plt.plot( U , mean_array , 'r', linewidth = 3.0 )
#b = plt.plot( U, lower_bound_scatter_band , 'k--', linewidth = 3.0 )
#c = plt.plot( U, upper_bound_scatter_band, 'k--' , linewidth = 3.0 )
plt.legend( ( z, a ), ( 'Realisationen', 'Erwartungswert' ) , 'upper right' )

plt.ylim( 0, 10 )
plt.show()

