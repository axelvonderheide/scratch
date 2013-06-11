
from enthought.traits.api import \
     Array, Bool, Callable, Enum, Float, HasTraits, Interface, implements, \
     Instance, Int, Trait, Str, Enum, Callable, List, TraitDict, Any, \
     on_trait_change, Tuple, WeakRef, Delegate, Property, cached_property, Dict, \
     Class, DelegatesTo

from enthought.traits.ui.api import \
     View, Item, Group
    
from numpy import \
     array, zeros, float_, dot, hstack, arange, argmin

from scipy.linalg import \
     det
     
from scipy.spatial.distance import \
     cdist     

from ibvpy.core.i_tstepper_eval import \
     ITStepperEval 
     
from ibvpy.core.tstepper_eval import \
     TStepperEval

from ibvpy.mats.mats_eval import \
    IMATSEval

from ibvpy.dots.dots_eval import \
    DOTSEval
    
from ibvpy.core.rtrace_eval import \
    RTraceEval

from ibvpy.fets.i_fets_eval import IFETSEval

#-------------------------------------------------------------------
# FETSEval - general implementation of the fe-numerical quadrature
#-------------------------------------------------------------------

class FETSEval( TStepperEval ):

    ngp_r = Int(0,label = 'Number of Gauss points in r-direction')

    n_gp = Property( depends_on = 'ngp_r')
    def _get_n_gp(self):
        return self.ngp_r * 10

class FETS2D4Q(FETSEval):

    # Integration parameters
    #
    ngp_r = Int(2)
    

class FETS2Drotsym( HasTraits ):

    parent_fets     = Instance( FETSEval )
    ngp_r           = DelegatesTo( 'parent_fets' )
        
if __name__ == '__main__':
    fe = FETS2Drotsym( fets = FETS2D4Q() )