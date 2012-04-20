'''
Created on Jul 14, 2009

@author: alexander
'''
from numpy import *
from scipy.linalg import *
from time import *


A2d = arange(4).reshape(2,2)
A3d = arange(9).reshape(3,3)

print 'A2d: \n', A2d
print 'A3d: \n', A3d

B = augment(A2d)
print 'B: \n', B