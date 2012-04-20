'''
Created on Aug 13, 2009

@author: rostislav
'''
from sympy import Symbol, cos
from sympy.mpmath import exp, e, pi
from symplot import mplot2d
from scipy import stats
from numpy import array, linspace
from matplotlib import pyplot as p

wb = Symbol('wb')
wl_arr = Symbol('wl_arr')
lR = Symbol('lR')
m = Symbol('m')
mb = Symbol('mb')
wsigma_b = Symbol('wsigma_b')
muT = Symbol('muT')
sigma_r = Symbol('sigma_r')

m = 5 * cos(mb)

mplot2d(m,(mb,0,6.1,100)) 








