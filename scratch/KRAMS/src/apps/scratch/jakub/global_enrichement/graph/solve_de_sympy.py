'''
Created on Apr 27, 2010

@author: jakub
'''
from sympy import *
x,y = symbols('xy')
f, g, h = map(Function, 'fgh')

pprint(f(x).diff(x, x) + f(x))

pprint(dsolve(f(x).diff(x, x) - f(x), f(x)))
