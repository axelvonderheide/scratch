'''
Created on 25.08.2011

@author: schmerl
'''
from enthought.traits.api import HasTraits, Property, cached_property, Event, \
    Array, Instance, Int, Directory, Range, on_trait_change, Bool, Trait, Constant, \
    Tuple, Interface, implements, Enum
from enthought.traits.trait_types import DelegatesTo
from enthought.traits.ui.api import Item, View, HGroup, RangeEditor
from numpy import loadtxt, min, array, arange, ones_like, cumsum, vstack, \
    hstack, sum, zeros_like, zeros, ones, where, unique, pi, invert, \
    prod,append
    
import math

class MyClass(HasTraits):
    arr = Array()
    #arr = (2,3),(4,5)
    arr
    c = arr.shape
    


if __name__ == '__main__':
    view1 = View(Item(name = 'arr'),
                 Item(name = 'c'))
    c = MyClass()
    c.configure_traits(view=view1)
    
        