

from dacwt import DAC

from enthought.traits.api import \
     Array, Bool, Callable, Enum, Float, HasTraits, Interface, implements, \
     Instance, Int, Trait, Str, Enum, Callable, List, TraitDict, Any, \
     on_trait_change, Tuple, WeakRef, Delegate, Property, cached_property, Dict

from enthought.traits.ui.api import \
     Item, View, HGroup, ListEditor, VGroup, Group

from enthought.traits.ui.menu import \
     NoButtons, OKButton, CancelButton, Action, CloseAction, Menu, \
     MenuBar, Separator

from math  import \
     pow, fabs

from numpy import \
     array, zeros, int_, float_, ix_, dot, linspace, hstack, vstack

from scipy.linalg import \
     inv, det

from ts_eval import \
     ITimeStepEval, TimeStepEval

from spatial_context import \
     ISDomain as ISD

from rv_eval import RTE

import time

class OnTrait(HasTraits):
    
    index1 = Int(-1)
    index2 = Int(-1)
    
    index3 = Property( Int, depends_on='index1,index2')
    def _set_index3(self,value):
        self._index3 = value
    @cached_property
    def _get_index3(self):
        return self.index1 * self.index2
    
    index4 = Property( Int, depends_on='index3' )
    @cached_property
    def _get_index4(self):
        return self.index3 * 2
    
ot = OnTrait( index1 = 3 )
ot.index2 = 8
print ot.index4
ot.index3 = 34
print ot.index4
ot.configure_traits()