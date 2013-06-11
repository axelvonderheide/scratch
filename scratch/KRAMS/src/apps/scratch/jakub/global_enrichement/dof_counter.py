'''
Created on Nov 24, 2009

@author: jakub
'''
from enthought.traits.api import \
     Array, Bool, Callable, Enum, Float, HasTraits, Interface, implements, \
     Instance, Int, Trait, Str, Enum, Callable, List, TraitDict, Any, Dict,\
     on_trait_change, Tuple, WeakRef, Delegate, Property, cached_property
     
class DofCounter(HasTraits):
    shape = Int(1)
    dof_on_edge = Int(2)
    
    dofs = Property(Int, depends_on = 'shape, shape')
    def _get_dofs(self):
        return self.shape*(1 + self.dof_on_edge-2)+1
    
if __name__ == '__main__':
    dc = DofCounter()
    dc.configure_traits()
        
