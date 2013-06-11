'''
Created on Apr 27, 2010

@author: jakub
'''
from sympy import *
x = symbols('x')

print integrate((cosh(x)/cosh(1.55))**2, (x, 0., 1.55))#.evalf(n=20)


