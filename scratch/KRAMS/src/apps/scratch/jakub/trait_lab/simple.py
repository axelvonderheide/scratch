'''
Created on Aug 3, 2009

@author: jakub
'''
from enthought.traits.api import \
    HasTraits, Instance, Property, Int, cached_property

class Simple(HasTraits):
    a = Int(1)
    
s = Simple()
s.configure_traits()