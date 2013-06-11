'''
Created on 24.09.2011

@author: axel
'''
from math import pi as Pi
import numpy as np
from scipy.stats import norm
from stats.pdistrib.sin2x_distr import sin2x
from scipy.stats import exponweib
from scipy.stats import weibull_min
from quaducom.resp_func.cb_short_fiber import CBShortFiber
from matplotlib import pyplot as plt

#set parameters
cbsf = CBShortFiber()
max = 0.01
nxz = 1000
U = np.linspace( 0, max, nxz )
anzahl_der_realisationen = 10
nx = 50
f = 1
length = 100
width = 30
height = 30
fiber_length = 9
area = width * height
r = 0.075
Vf = 3
def H( x ):
    return np.sign( np.sign( x ) + 1. )

start = 0.0001
end = 0.99999
dp = 3000
############################################
Le = fiber_length / 2
specimen_volume = length * width * height
no_of_fibers_in_specimen = ( specimen_volume * Vf / 100 ) / ( Pi * r ** 2 * fiber_length )
prob_crackbridging_fiber = .5 * fiber_length / length
mean_nof = prob_crackbridging_fiber * no_of_fibers_in_specimen
var_nof = ( prob_crackbridging_fiber * no_of_fibers_in_specimen * ( 1 - prob_crackbridging_fiber ) )
stdev_nof = var_nof ** 0.5



x_array = np.linspace( start, end, dp )

def pdf_tau( x ) :
    return  weibull_min( 8.50666480278, 2.49019483074 ).pdf( x )


def ppf_tau ( x ):
    return weibull_min( 8.50666480278, 2.49019483074 ).ppf( x )

def pdf_phi ( x ):
    return sin2x.pdf( x )

def ppf_phi ( x ):
    return sin2x.ppf( x )

'''
def pdf_f ( x ):
    return  norm( 0.87, 0.2 ).pdf( x )

def ppf_f ( x ):
    return  norm( 0.87, 0.2 ).ppf( x )
'''

def pdf_le ( x ):
    return 1. / Le * ( x ** 0 )



def ppf_le ( x ):
    return x * Le

def pdf_n( x ):
    return  norm( mean_nof, stdev_nof ).pdf( x )

def ppf_n ( x ):
    return norm( mean_nof, stdev_nof ).ppf( x )


#Realisationen der Zufallsparameter

real_tau_or = np.linspace( ppf_tau( start ), ppf_tau( end ), nx )
print real_tau_or
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
    #step_three = np.trapz( step_two, x_phi )
    step_four = np.trapz( step_two, x_le )
    mean_matrix_list.append( step_four )

    stdev_to_be_integrated = ( pullouts - step_four ) ** 2 * pdf_matrix
    step_one1 = np.trapz( stdev_to_be_integrated, x_tau )
    step_two2 = np.trapz( step_one1, x_phi )
    #step_three3 = np.trapz( step_two2, x_phi )
    step_four4 = np.trapz( step_two2, x_le )
    variance_list.append( step_four4 )


mean_array = np.array( mean_matrix_list ) * mean_nof
var_fiber_array = np.array( variance_list )
var_array = var_fiber_array * mean_nof + stdev_nof * mean_array
stdev_array = var_array ** 0.5
lower_bound_scatter_band_pre = mean_array - stdev_array
lower_bound_scatter_band = H( lower_bound_scatter_band_pre ) * lower_bound_scatter_band_pre
upper_bound_scatter_band = mean_array + stdev_array
lower_twice_stdev = mean_array - stdev_array * 2
upper_twice_stdev = mean_array + stdev_array * 2
lower_triple_stdev = mean_array - stdev_array * 3
upper_triple_stdev = mean_array + stdev_array * 3
#print np.max( mean_array )
print stdev_array[np.argmax( mean_array )]
farbe = 'black'
#plt.plot( U, lower_twice_stdev, 'k--', U, upper_twice_stdev, 'k--' )
#plt.plot( U, lower_triple_stdev, 'k--', U, upper_triple_stdev, 'k--' )
plt.title( 'Kraft-Rissoeffnungskurve', fontsize = 16 )
plt.xlabel( 'Rissoeffnung $w$ in [mm]', fontsize = 16 )
plt.ylabel( '$\sigma_f$[N/mm$^2$]', fontsize = 16 )
one_stdev = plt.fill_between( U, lower_triple_stdev / area, upper_triple_stdev / area  , color = farbe, alpha = 0.4 )
twice_stdev = plt.fill_between( U, lower_twice_stdev / area, upper_twice_stdev / area , color = farbe, alpha = 0.4 )
triple_stdev = plt.fill_between( U, lower_bound_scatter_band / area, upper_bound_scatter_band / area , color = farbe, alpha = 0.4 )
plt.legend( ( one_stdev, twice_stdev, triple_stdev ), ( '99 Umgebung', '95 Umgebung', '68 Umgebung' ), 'lower right' )

w = np.linspace( 0, max, nxz )
U = U.reshape( nxz, 1 )
for i in range( anzahl_der_realisationen ):
    nof = round( ppf_n( np.random.rand( 1 ) ) )
    tau_realisationen = ppf_tau( np.random.rand( nof ) ).reshape( 1, nof )
    phi_realisationen = ppf_phi( np.random.rand( nof ) ).reshape( 1, nof )
    f_realisationen = ppf_phi( np.random.rand( nof ) ).reshape( 1, nof )
    phi_realisationen = ppf_phi( np.random.rand( nof ) ).reshape( 1, nof )
    le_realisationen = ppf_le( np.random.rand( nof ) ).reshape( 1, nof )
    pullouts = cbsf( U, tau_realisationen, 0, 0.15, 200000, le_realisationen, phi_realisationen, f, 0, 0, 1e15 )
    results = np.sum( pullouts , axis = 1 )
    plt.plot( w, results / area, 'b' )
bum = plt.plot( w, results / area , 'b' )
mean_plot = plt.plot( U , mean_array / area , 'r', linewidth = 3.0 )
#b = plt.plot( U, lower_bound_scatter_band / area, 'k--', linewidth = 3.0 )
#c = plt.plot( U, upper_bound_scatter_band / area, 'k--' , linewidth = 3.0 )
flu = [-1, -2]
flu1 = [3, 4]
z = plt.fill( flu, flu1, color = farbe, alpha = 1 )
a = plt.fill( flu, flu1, color = farbe, alpha = .8 )
b = plt.fill( flu, flu1, color = farbe, alpha = .4 )
plt.legend( ( bum, mean_plot, z, a, b ), ( 'Realisationen', '$E_\sigma(\sigma)$', '$\sigma$' , '2$\sigma$', '3$\sigma$' ), 'lower right' )

plt.ylim( 0 )
plt.xlim ( 0, max )
plt.show()

