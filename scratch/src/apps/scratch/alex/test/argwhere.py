'''
Created on Jul 2, 2009

@author: alexander
'''
from numpy import *

x = array([13,14,15,165,16,16])
print 'x',x

b = array([1,2,3,4,16,16])
print 'indexing', x[[0,2]]

ca = x.all()
print 'ca', ca

c = b == x 
print 'where', c.all()

#
#y = array([3,4,5,65,6,6])
#print 'y', y
#
#b = x[where(x<=30.)[0]]
#print 'b',b
#
#idx_emax = len(b)
#c = y[ :idx_emax ]
#print 'c',c
#
#
#
#b=argwhere(x<=1.)
#print 'b',b[0]
#
#print 'x[b[:]]', x[b[:]]
