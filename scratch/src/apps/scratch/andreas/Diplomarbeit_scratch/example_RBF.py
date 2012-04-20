'''
Created on Sep 15, 2009

@author: Andreas
'''
from enthought.traits.api import \
    HasTraits, Float, Array, implements, Property, cached_property, Instance,\
    Int

from mathkit.mfn.mfn_line.mfn_line import MFnLineArray

from numpy import \
    array, tensordot, dot,zeros, c_, ix_, mgrid, arange

from math import sqrt, asin, acos

# Interpolation
from scipy.interpolate import Rbf

from scipy import rand

from simiter.sim_pstudy import ISimModel, SimOut, SimPStudy
 
 
 
def test(points):
    
    # create some points
    
    x = rand(10)
            
    y = rand(10)
    
    z = rand(10)
    
    d = x**2+y**2-1
    
    rbf = Rbf(x, y, z )#function = 'gaussian')
    
    print "done"

    return
test(10)