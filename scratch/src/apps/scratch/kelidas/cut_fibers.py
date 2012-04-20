'''
Created on 24.3.2011

@author: Kelidas
'''


from sympy import sin, Symbol, cos, pprint, integrate, Rational, N, evalf
from mpmath import pi
ell = Symbol( 'ell' )
phi = Symbol( 'phi' )
L = Symbol( 'L' )
x = Symbol( 'x' )
le = Symbol( 'le' )

fphi = sin( phi )
fx = 1 / L
ff = fphi * fx # joint density function
p = Rational( 1, 2 ) * ell / L # probability 

#pprint( ff )

# indicator function change integration domain
res = integrate( ff , ( x, -Rational( 1, 2 ) * ell * cos( phi ),
                         Rational( 1, 2 ) * ell * cos( phi ) ) ) / p


pprint( res )
