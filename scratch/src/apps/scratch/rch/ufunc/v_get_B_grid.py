
from enthought.traits.api import \
     Array, Bool, Callable, Enum, Float, HasTraits, Interface, implements, \
     Instance, Int, Trait, Str, Enum, Callable, List, TraitDict, Any, \
     on_trait_change, Tuple, WeakRef, Delegate, Property, cached_property

from ibvpy.fets.fets2D.fets2D4q import FETS2D4Q
from ibvpy.fets.fets2D.fets2D4q16u import FETS2D4Q16U
from numpy import frompyfunc, apply_along_axis, \
    apply_over_axes, sum, array, expand_dims, copy, zeros, vectorize, \
    dot, tensordot, c_

from scipy.linalg import det, inv


#------------------------------------------------------------------------------
# Run a loop over the elements and put the data into the 
#------------------------------------------------------------------------------
def apply_on_ip_grid( fn, X_el, ip_array, ip_state = None ):
    '''
    Apply the function fn over the first dimension of the array.
    
    '''
    out_single = fn( ip_array[0], X_el[0] )
    out_grid_shape = ( X_el.shape[0], ip_array.shape[0], ) + out_single.shape
    out_grid = zeros( out_grid_shape )

    for el in range( X_el.shape[0] ):
        for ip in range( ip_array.shape[0] ):
            out_grid[ el, ip, ... ] = fn( ip_array[ip], X_el[el] )

    return out_grid


def get_B_mtx_grid( domain, fets_eval ):
    # elementwise specification of the coordinates
    X_el = domain.elem_X_map  # X[ el, el_node, coord ]
    ip_array = fets_eval.ip_array

    B_mtx_grid = apply_on_ip_grid( fets_eval.get_B_mtx, X_el, ip_array )
    return B_mtx_grid

if __name__ == '__main__':
    # Get the global coordinates of all integration points in the domain.
    #
    from ibvpy.mesh.fe_grid import FEGrid
    # Discretization
    domain = FEGrid( coord_min = ( 0., -1., 0. ),
                           coord_max = ( 1., -0., 0. ),
                           shape = ( 2, 1 ),
                           n_nodal_dofs = 2,
                           dof_r = [[-1, -1], [1, -1], [1, 1], [-1, 1]],
                           geo_r = [[-1, -1], [1, -1], [1, 1], [-1, 1]] )

    fets_eval = FETS2D4Q()
    B_mtx_grid = get_B_mtx_grid( domain, fets_eval )
    print 'B_mtx_grid FETS2D4Q'
    print B_mtx_grid

    fets_eval = FETS2D4Q16U()
    B_mtx_grid = get_B_mtx_grid( domain, fets_eval )
    print 'B_mtx_grid FETS2D4Q'
    print B_mtx_grid
