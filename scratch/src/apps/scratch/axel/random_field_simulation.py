'''
Created on 24.06.2011

@author: axel
'''
from numpy import linspace
from math import e
from scipy.linalg import toeplitz
from numpy.linalg import eig, inv
from numpy import dot, transpose, ones, array, sum, var , eye , real, mean
from scipy.stats import norm
from numpy.random import shuffle
from matplotlib import pyplot as p

def acor( dx ):
    return e ** ( -( dx / lcorr ) ** 2 )

lcorr = 0.1
L = 1000.
nd = 1000.
xgrid = linspace( 0, L, nd )
Rdist = toeplitz( xgrid, xgrid )
R = acor ( Rdist )

_lambda, phi = eig( R )

#print _lambda , phi

nsim = 1.
randsim = linspace( 0, 1, nd + 1 ) - 0.5 / ( nd )
randsim = randsim[1:]
shuffle( randsim )
xi = transpose( ones( ( nsim, nd ) ) * array( [ norm().ppf( randsim ) ] ) )
#print phi.shape, _lambda.shape, xi.shape
lam = eye( nd ) * _lambda
H = dot( dot( phi, ( lam ) ** 0.5 ), xi )
H = H.real

p.plot( xgrid, H )
p.show()

