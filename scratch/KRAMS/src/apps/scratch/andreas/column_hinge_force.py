
'''
Created on Jul 30, 2010

@author: abach
'''
from enthought.traits.api import \
    HasTraits, Directory, List, Int, Float, Any, Enum, \
    on_trait_change, File, Constant, Instance, Trait, \
    Array, Str, Property, cached_property, WeakRef, \
    Dict, Button, Color, Bool

from scipy import integrate
from numpy import array, linspace, ones_like

from math import pi

steps = int ( 1000 )
x_i = linspace( 0, 3, steps ) # m
M = linspace( 1, 4, steps )  #m
I_i = ( ( -0.1 * x_i / 3.0 + 0.45 ) ** 4 ) / 12  # m^4
E = 38800  # MN/m^2

Func = M ** 2 / ( E * I_i ) # 1/ m
Integral = integrate.trapz( Func, x_i )
Force = 0.0052 / Integral
print Force

