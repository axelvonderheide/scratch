''' 
The first example evaluates the term in the brackets first and because 1 (=numerator of the fraction) 
and ndiv are  of type integer 1/ndiv evaluates to 0!

In the second example the type of the result is derived from 0.7 
so that the result is also of type float!
'''
from numpy import *
from enthought.traits.api import Array, Bool, Callable, Enum, Float, Event, HasTraits, \
                                 Instance, Int, Trait, ToolbarButton, Button, on_trait_change, \
                                 Property, cached_property
                                 
#import matplotlib.pyplot as plt
#plt.plot([1,2,3],[2,3,4])
#plt.ylabel('XXX')
#plt.show()


a = 3*('XXX'+',')
print 3*a
#                                 
#                                 
#                                 
#e_max = 0.030303030303
#Epp = 0.0
#
#print 'sqrt', sqrt( Epp / e_max * exp( -(e_max-Epp)/Efp ))
#
#
#x = arange(10)
#print 'x', x
#
#
#
#a = x[:3]
#b = x[4:]
#print 'a', a
#print 'b', b
#
#c = hstack([a,b])
#print 'c', c
#
#d = delete(x,-2, None)
#print 'd', d
#
#a = array([0.,1.])
#print 'a', a
#
#b = (x == a).all() 
#print 'b', b
#
#
#
#
#x_.all() == array([0., 1.]).all():
#
#
#
#xdata = array([0.0, 1.0, 2.0])
#ydata = array([0.0, 4.0, 3.0])
#print 'xdata: ', xdata
#print 'ydata: ', ydata
#
#x = 0.2
#print 'x: ', x
#
#x2idx = xdata.searchsorted(x)
#print 'x2idx: ', x2idx
#
#if x2idx == len(xdata):
#    x2idx -= 1
#
#x1idx = x2idx - 1
#print 'x1idx: ', x1idx
#
#x1 = xdata[ x1idx ]
#print 'x1: ', x1
#
#x2 = xdata[ x2idx ]
#print 'x2: ', x2
#
#dx = x2 - x1
#print 'dx: ', dx
#
#y1 = ydata[ x1idx ]
#print 'y1: ', y1
#
#y2 = ydata[ x2idx ]
#print 'y2: ', y2
#
#dy = y2 - y1
#
#y = y1 + dy / dx * (x - x1)
#
#print 'y: ', y
#
#
#
#
#
#
#
#
#
#
#
#
#x = lambda angle: 1., 1, 1
#print x(1)
#
#A = arange(10)
#B = array([10,0,16,66,7,8,2,3,0,0])
#print 'A_old', A
#print 'B', B
#
#bool_arr = B >= A
#A[bool_arr] = B[bool_arr]
#print 'bool_arr', bool_arr
#print 'A_new', A
#
#
#
##ndiv = 4
##for i in range(4):
##    print 0.7 * (1/ndiv) * i
##    
##print '\n'
##
##ndiv = 4
##for i in range(4):
##    print 0.7 * 1/ndiv * i
##    
##print 'range(4)[1]:', type(range(4)[1])
##
##a = 0.7 * (1/4)
##b = 0.7 * 1/4 
##
##print a 
##print b 