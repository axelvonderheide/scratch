'''
Created on Apr 6, 2009

@author: jakub
'''
from enthought.traits.api import \
    Array, Bool, Callable, Enum, Float, HasTraits, \
    Instance, Int, Trait, Str, Enum, Callable, List, TraitDict, Any, \
    Dict, Property, cached_property, WeakRef, Delegate, Button, \
    Constant, HasTraits, Tuple

from ibvpy.rtrace.rt_domain import RTraceDomain
from ibvpy.rtrace.rt_domain_list import RTraceDomainList

class RTraceDL(HasTraits):
    label = Property
    def _get_label(self):
        return (1,2,3)

class RTraceMy(RTraceDL):
    
    def show(self):
        l = super(RTraceMy,self)._get_label() 
        #l = self.label
        print l

rt = RTraceMy()
rt.show()