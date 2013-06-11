
from enthought.traits.api import \
     Array, Bool, Callable, Enum, Float, HasTraits, Interface, implements, \
     Instance, Int, Trait, Str, Enum, Callable, List, TraitDict, Any, \
     on_trait_change, Tuple, WeakRef, Delegate, Property, cached_property, Dict

from enthought.traits.ui.api import \
     Item, View

from enthought.traits.ui.menu import \
     OKButton, CancelButton
     

from numpy import \
     zeros, float_, ix_

from core.ts_eval import \
     ITStepperEval, TStepperEval

from core.rtrace_eval import RTraceEval
from core.fets_eval import IFETSEval

from core.dots_eval import DOTSEval

from scipy import sparse

def flatten(x):
    """flatten(sequence) -> list

    Returns a single, flat list which contains all elements retrieved
    from the sequence and all recursively contained sub-sequences
    (iterables).

    Examples:
    >>> [1, 2, [3,4], (5,6)]
    [1, 2, [3, 4], (5, 6)]
    >>> flatten([[[1,2,3], (42,None)], [4,5], [6], 7, MyVector(8,9,10)])
    [1, 2, 3, 42, None, 4, 5, 6, 7, 8, 9, 10]"""

    result = []
    for el in x:
        #if isinstance(el, (list, tuple)):
        if hasattr(el, "__iter__") and not isinstance(el, basestring):
            result.extend(flatten(el))
        else:
            result.append(el)
    return result


#-----------------------------------------------------------------------------
# Integrator for a simple 1D domain.
#-----------------------------------------------------------------------------

class SpDOTSEval( DOTSEval ):
    '''
    Domain with uniform FE-time-step-eval.
    '''
    implements( ITStepperEval )

    fets_eval = Instance( IFETSEval )
    some = Str('')

    def new_cntl_var(self):
        return zeros( self.domain.n_dofs, float_ )
    
    def new_resp_var(self):
        return zeros( self.domain.n_dofs, float_ )
    
    def get_state_array_size( self ):

        # The overall size is just a n_elem times the size of a single element
        #
        shape = sctx.sdomain.shape
        #self.e_arr_size = shape * self.fets_eval.get_state_array_size( )
        self.e_arr_size = self.fets_eval.get_state_array_size( )
        dots_arr_size = shape * self.e_arr_size
        return dots_arr_size

    def setup( self, sctx ):
        self.domain = sctx.sdomain
        
        ndofs = self.domain.n_dofs
        self.K = sparse.lil_matrix((ndofs, ndofs), float_ )
        self.K_temp = sparse.lil_matrix((ndofs, ndofs), float_ )
        self.F_int = zeros( ndofs, float_ )

        e_arr_size      = self.e_arr_size

        # Run the setup of sub-evaluator
        #
        for elem in sctx.sdomain.elements:
            sctx.elem = elem
            e_id = elem.id_number
            sctx.elem_state_array = sctx.state_array[ e_id * e_arr_size : (e_id + 1) * e_arr_size ]
            self.fets_eval.setup( sctx )
    
    def get_corr_pred( self, sctx, u, tn, tn1 ):
        ndofs = self.domain.n_dofs
        #self.K.data[::] = 0.0
        self.K = sparse.lil_matrix((ndofs, ndofs), float_ )
        self.F_int[:]    = 0.0
        e_arr_size      = self.e_arr_size

        for elem in sctx.sdomain.elements:
            e_id = elem.id_number
            ix = elem.get_dof_map()
            sctx.elem = elem
            sctx.elem_state_array = sctx.state_array[ e_id*e_arr_size : (e_id+1)*e_arr_size ]
            sctx.X = elem.get_X_mtx()
            f, k = self.fets_eval.get_corr_pred( sctx, u[ix_(ix)], tn, tn1 )
            #self.K_temp.data[:][:] = 0.
            self.K_temp = sparse.lil_matrix((ndofs, ndofs), float_ )
            a = 0
            for i in ix:
                self.K_temp.rows[i] = ix
                self.K_temp.data[i][:] = k[a][:]
                a =+1
                #print K_temp
            self.K = self.K + self.K_temp
            self.F_int[ ix_(ix) ] += f


        return self.F_int, self.K

    def map_u(self, sctx, U):
        ix = sctx.elem.get_dof_map()
#        sctx.r = fets_eval.map_to_local( sctx.elem, sctx.X )
        u = U[ix]
        return u
    
    rte_dict = Property( Dict, depends_on = 'fets_eval')
    @cached_property
    def _get_rte_dict(self):
        rte_dict = {}
        for key, eval in self.fets_eval.rte_dict.items():
            rte_dict[ key ] = RTraceEvalUDomainFieldVar( name = key, 
                                                  u_mapping = self.map_u, 
                                                  eval = eval,
                                                  fets_eval = self.fets_eval )
        return rte_dict
      
    traits_view = View( Item('fets_eval', style = 'custom', show_label = True ),
                        Item('some', style = 'custom', show_label = False ),
                        resizable = True,
                        height = 0.8,
                        width = 0.8,
                        buttons = [OKButton,CancelButton],
                        kind = 'subpanel',
                        scrollable = True,
                        )



class RTraceEvalUDomainFieldVar(RTraceEval):
    fets_eval = WeakRef( IFETSEval )
    
    # @TODO Return the parametric coordinates of the element covering the element domain
    #
    vtk_r_arr = Delegate('fets_eval')
    field_entity_type = Delegate('fets_eval')
    dim_slice = Delegate('fets_eval')
    get_vtk_r_glb_arr = Delegate('fets_eval')
    n_vtk_r = Delegate('fets_eval')
    field_faces = Delegate('fets_eval')
    field_lines = Delegate('fets_eval')
    get_state_array_size = Delegate('fets_eval')
    
    # @TODO Return the entity lists of the element used for element visualization
    #
    def _get_points(self, X):
        return fets_eval.get_points()
    
    
