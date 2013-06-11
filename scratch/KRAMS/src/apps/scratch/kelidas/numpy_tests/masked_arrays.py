'''
Created on Jan 25, 2011

@author: kelidas
'''

from numpy import *
import numpy.ma as ma




a = array( [[1, 2, 3], [4, 5, 6], [1, 5, 3]] )
m1 = array( [[True, False, False], [False, True, False], [False, False, False]], dtype = bool )
m2 = array( [[False, False, True], [False, True, False], [False, False, False]], dtype = bool )

# creating masked array
# value at position 'True' will be hidden
a = ma.masked_array( a, mask = m1, fill_value = -100 )
b = ma.masked_array( a, mask = m2, fill_value = -100 )


print 'masked array a', a
print 'masked array b', b

print 'array without mask -- current data', a.data
print 'array mask', a.mask
print 'array (1D) of non-masked values', a.compressed()
print 'unmasked array, masked positions are filled with fill_value', a.filled()



# evaluate only values non-masked in both arrays, other values in result array will be masked
print 'sum of two arrays with different masks', a + b
print 'product of two arrays with different masks', a * b

# in numpy.ma there is many functions similar to pure numpy package
# numpy.mean and ma.mean give the same result, but this needn't be satisfied in all functions 
print 'mean value in the row', ma.mean( a, axis = 1 )

print 'sum values in the row', ma.sum( a, axis = 1 )








