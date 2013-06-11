'''
Created on 30.09.2011

@author: axel
'''

from math import pi
import numpy as np

x = np.linspace( -0.5, 0.5, 50 )
x = x.reshape( 1, 50 )
x_integ = np.linspace( -0.5, 1, 50 )
phi = np.linspace( 0, pi / 2, 50 )
phi = phi.reshape( 50, 1 )
phi_integ = np.linspace( 0, pi / 2, 50 )
func = ( 0.5 * np.cos( phi ) - abs( x ) )
mean = np.trapz( func, phi_integ )
mean1 = np.trapz( mean, x )
print mean1
