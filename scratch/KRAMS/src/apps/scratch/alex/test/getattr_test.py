'''
Created on Jul 16, 2010

@author: alexander
'''
from enthought.traits.api import \
    HasTraits, Directory, List, Int, Float, Any, Enum, \
    on_trait_change, File, Constant, Instance, Trait, \
    Array, Str, Property, cached_property, WeakRef, \
    Dict, Button, Color, Bool, Trait

from numpy import *

#from itertools import \
#    permutation, combination


def _product( args ):
    """
    Get all possible permutations of the security factors
    without changing the order of the loading cases.
    The method corresponds to the build-in function 'itertools.product'.
    Instead of returning a generator object a list of all 
    possible permutations is returned. As argument a list of list 
    needs to be defined. In the original version of 'itertools.product' 
    the function takes a tuple as argument ("*args"). 
    """
    pools = map( tuple, args ) #within original version args defined as *args
    result = [[]]
    for pool in pools:
        result = [x + [y] for x in result for y in pool]
    return result


print _product( [range( 3 ), range( 3 )] )
sl = ogrid[ 0:3, 0:3 ]
a = arange( 9 ).reshape( 3, 3 )
print a[sl]


