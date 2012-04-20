#-------------------------------------------------------------------------------
#
# Copyright (c) 2009, IMB, RWTH Aachen.
# All rights reserved.
#
# This software is provided without warranty under the terms of the BSD
# license included in simvisage/LICENSE.txt and may be redistributed only
# under the conditions described in the aforementioned license.  The license
# is also available online at http://www.simvisage.com/licenses/BSD.txt
#
# Thanks for using Simvisage open source!
#
# Created on Apr 12, 2010 by: rch

# define a three-dimensional array

from numpy import arange, where, broadcast_arrays, ones_like, ones
from scipy.integrate import trapz, simps, romb

arr1 = arange( 2*3*4, dtype = 'float' ).reshape((2,3,4))
arr2 = arange( 3*4, dtype = 'float' ).reshape((3,4))


print 'arr1\n', arr1
print 'arr2\n', arr2

arr1_x_arr2 = arr1 * arr2
print 'arr1 * arr2\n', arr1_x_arr2

i0_arr1_x_arr2 = trapz( trapz( trapz( arr1_x_arr2, axis = 0 ), axis = 0), axis = 0 )
i1_arr1_x_arr2 = trapz( trapz( trapz( arr1_x_arr2, axis = 2 ), axis = 1), axis = 0 )
print 'int( arr1 * arr2 )', i0_arr1_x_arr2, i1_arr1_x_arr2

a1 = arange( 3, dtype = float ).reshape((1,3))
a2 = arange( 3, dtype = float ).reshape((3,1))
print 'a1', a1
print 'a2', a2
a1_ = a1 * ones_like( a1 * a2 )
print a1_
a1_[ a1 <= a2 ] = 100
print a1_
 
