'''
Created on 12.06.2013

@author: acki
'''
import numpy as np
from matplotlib import pyplot as plt
from stats.spirrid.rv import RV

a = np.array( [1, 2, 3] ).reshape( 3, 1 )
b = np.array( [4, 5, 6] ).reshape( 1, 3 )
c = np.array( [2, 3, 100] ).reshape( 3, 1 )
matrix = a * b

print a * b, a * b + c, np.sum( a * b, 0 )
