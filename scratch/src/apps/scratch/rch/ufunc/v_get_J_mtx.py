
from enthought.traits.api import \
     Array, Bool, Callable, Enum, Float, HasTraits, Interface, implements, \
     Instance, Int, Trait, Str, Enum, Callable, List, TraitDict, Any, \
     on_trait_change, Tuple, WeakRef, Delegate, Property, cached_property

from ibvpy.fets.fets2D.fets2D4q import FETS2D4Q
from numpy import frompyfunc, apply_along_axis, \
    apply_over_axes, sum, array, expand_dims, copy, zeros, vectorize, \
    dot, tensordot, c_

from scipy.linalg import det, inv

fets_eval = FETS2D4Q()

#------------------------------------------------------------------------------
# Run a loop over the elements and put the data into the 
#------------------------------------------------------------------------------
def apply_along_first_axis( fn, inp ):
    '''
    Apply the function fn over the first dimension of the array.
    
    In 
    '''
    # Get the shape of the output by taking the first dimension
    # of the input with the shape of the result of fn
    #
    out_shape = ( inp.shape[0], ) + fn( inp[0] ).shape

    # setup the output as zero valued array
    #
    out = zeros( out_shape, dtype = 'float_' )

    # run a loop over the first dimension and store the result 
    # of fn in the array
    #
    for i, inp_entry in enumerate( inp ):
        out[i, ...] = fn( inp_entry )
    return out

def get_N_geo_mtx( r_pnt ):
    '''Reduce the dimension of the N matrix (actually, it should be
    reduced directly in the implementation as it can be always 
    considered a vector.
    '''
    N_geo_mtx = fets_eval.get_N_geo_mtx( r_pnt )
    return N_geo_mtx.reshape( ( N_geo_mtx.shape[1], ) )

# Evaluate the value of N_r at all integration points 
# of a single element
#
N_ip = apply_along_first_axis( get_N_geo_mtx, fets_eval.ip_array )
print 'N_ip'
print N_ip  # [ ip_node, el_node ]  # el_node is the element corner node

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

# elementwise specification of the coordinates
X_el = domain.elem_X_map  # X[ el, el_node, coord ]

print 'X_ip'
print X_el
print X_el.shape

# Be careful - the index ordering gets change in the following way
#
# X_ip[ el, coord, ip_node ] = X_mtx[ el, el_node, coord ] * N_ip[ ip_node, el_node ]
# That's the reason why to swap the axes after the call
#
X_ip = tensordot( X_el, N_ip, axes = ( -2, -1 ) ).swapaxes( 1, 2 )
print 'X_ip'
print X_ip

# An alternative would be
#
# X_ip[ ip_node, el, coord ] = N_ip[ ip_node, el_node ] * X_el[ el, el_node, coord ]
# Thus, here axes 0,1 must be swapped

X_ip = tensordot( N_ip, X_el, axes = ( -1, -2 ) ).swapaxes( 0, 1 )
print 'X_ip'
print X_ip

# The above call is equivalent to the call
#
X_ip = dot( N_ip, X_el ).swapaxes( 0, 1 )
print 'X_ip'
print X_ip

# get the list of points

pnts = c_[ X_ip[..., 0].flatten(), X_ip[..., 0].flatten()]
print 'pnts'
print pnts

###### Now evaluate the jacobi matrix for each ip

dN_ip = apply_along_first_axis( fets_eval.get_dNr_geo_mtx, fets_eval.ip_array )
print 'dN_ip[ ip_node, coord, el_node ]'
print dN_ip

J_ip = tensordot( X_el, dN_ip, axes = ( -2, -1 ) ).swapaxes( 1, 2 )
print 'J_ip[ el, coord, ip_node, coord ] = X_el[ el, el_node, coord ] * dN_ip[ ip_node, coord, el_node ]'
#print J_ip

J_ip = dot( dN_ip, X_el ).swapaxes( 1, 2 ).swapaxes( 0, 1 )
print 'J_ip[ ip_node, coord, el, coord ] = dN_ip[ ip_node, coord, el_node ] * X_el[ el, el_node, coord ]'
#print J_ip

for el in range( J_ip.shape[0] ):
    for ip in range( J_ip.shape[1] ):
        print det( J_ip[ el, ip, :, :] )
        print inv( J_ip[ el, ip, :, :] )
