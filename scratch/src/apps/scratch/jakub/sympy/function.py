'''
Created on Aug 10, 2009

@author: jakub
'''
from sympy import symbols, Matrix
from sympy.functions import sign
from sympy.functions.special.delta_functions import Heaviside, DiracDelta



x,x_i,Xi = symbols('x,x_i,Xi') 

def Jump(x,x_i,Xi):
    return 1./2.*(sign(x-Xi)-sign(x_i-Xi))


print Jump(x,0., 0.5)