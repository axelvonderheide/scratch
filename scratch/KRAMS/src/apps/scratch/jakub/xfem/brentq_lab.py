from scipy.optimize import brentq

from numpy import \
     array, zeros, int_, float_, ix_, dot, linspace, hstack, vstack, arange, \
     identity, unique, average, frompyfunc, linalg, sign

from enthought.traits.api import \
     Array, Bool, Callable, Enum, Float, HasTraits, Interface, implements, \
     Instance, Int, Trait, Str, Enum, Callable, List, TraitDict, Any, \
     on_trait_change, Tuple, WeakRef, Delegate, Property, cached_property


def ls_function( r, s): 
    return r -0.2#Y-0.2*X#temp

def ls_fn_r(r,s):
    return ls_function(r,s)

def ls_fn_s(s,r):
    return ls_function(r,s)

pt_sets = Property(List)
def _get_pt_sets():
    inter_pts = []
    for c_coord in [-1,1]:
        args = (c_coord)
        s_coord = _get_intersect_pt(ls_fn_s, args)
        r_coord = _get_intersect_pt(ls_fn_r, args)
        if s_coord:
            inter_pts.append([c_coord,s_coord])
        if r_coord:
            inter_pts.append([r_coord,c_coord])
    return inter_pts

def _get_intersect_pt(fn, args):
    try:
        return brentq(fn,-1,1,args=args)
    
    except ValueError:
        return

def _get_ls_nodes():
    coords = array([[-1,-1],[1,-1],[1,1],[-1,1]])
    r,s = coords.T
    ls_fn = frompyfunc( ls_function, 2, 1 )
    ls_vals = ls_fn(r,s)
    pos_arr = coords[ls_vals[:] > 0.]
    neg_arr = coords[ls_vals[:] < 0.]
    return pos_arr, neg_arr
    


    
print _get_pt_sets()
print _get_ls_nodes()