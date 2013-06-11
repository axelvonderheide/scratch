'''
Created on Apr 3, 2012

@author: rostar
'''
from matplotlib import pyplot as plt
from scipy.stats import weibull_min, uniform
import numpy as np
from stats.misc.random_field.random_field_1D import RandomField
from stats.spirrid import make_ogrid
#from enthought.mayavi import mlab as m
from math import pi
from stats.pdistrib.sin2x_distr import sin2x
            
l = np.linspace(0.0, 1., 100)
phi = np.linspace(0.0,pi/2, 100)

# joint probability density function of embedded lengths and inclination angles
#e = make_ogrid([l, phi])
#m.surf(e[0], e[1], l*np.cos(phi[:,np.newaxis]))
#m.surf(e[0], e[1], 1./2*np.ones_like(l)*np.sin(2.*phi[:,np.newaxis]))
#m.show()

# exact solution of number of fibers with distance from crack
sim = 10000.

def mc(z, uhel, delka):
    return np.sum(delka * np.cos(uhel) > z)/sim

n_list = []
z_arr = np.linspace(0,1,50)
for z in z_arr:
    n_list.append(mc(z,sin2x.ppf(np.random.rand(sim)),
                     uniform(loc = 0, scale = 1).ppf(np.random.rand(sim))))

def exact(lf, z):
    return (lf-2*z)**2/lf**2

def power(z):
    return np.ones_like(z) - 2*z + z**2

plt.plot(z_arr,n_list, color = 'red', lw = 3, ls = 'dashed', label = 'MC')
plt.plot(z_arr,exact(2.,z_arr), color = 'black', lw = 1, label = 'analytical')
plt.plot(z_arr,power(z_arr), color = 'green', lw = 2, label = 'power series')
plt.legend(loc = 'best')
plt.ylim(0,1.1)
plt.show()