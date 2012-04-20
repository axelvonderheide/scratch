
from enthought.traits.api import \
     Array, Bool, Callable, Enum, Float, HasTraits, Interface, implements, \
     Instance, Int, Trait, Str, Enum, Callable, List, TraitDict, Any, \
     on_trait_change, Tuple, WeakRef, Delegate, Property, cached_property

from enthought.traits.ui.api import \
     Item, View, HGroup, ListEditor, VGroup, Group

from enthought.traits.ui.menu import \
     NoButtons, OKButton, CancelButton, Action, CloseAction, Menu, \
     MenuBar, Separator

from numpy import \
     array, zeros, int_, float_, ix_, dot, linspace, hstack, vstack, arange, \
     identity, unique, average, frompyfunc, linalg, sign

from scipy.linalg import \
     inv
     
from scipy.optimize import brentq

from fets_eval_enr import FETSEvalEnr
from ibvpy.fets.fets2D.fets2D4q import FETS2D4Q
from ibvpy.fets.fets2D.fets2D4q8u import FETS2D4Q8U
from ibvpy.fets.fets2D.fets2D4q9u import FETS2D4Q9U
from ibvpy.fets.fets2D.fets2D9q import FETS2D9Q

from XFETSEval import XFETSEval

#-----------------------------------------------------------------------------------
# FETS2D4Q - 4 nodes isoparametric quadrilateral element (2D, linear, Lagrange family)    
#-----------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------
# Element Information: 
#-----------------------------------------------------------------------------------
#
# Here an isoparametric element formulation is applied.
# The implemented shape functions are derived based on the 
# ordering of the nodes of the parent element defined in 
# '_node_coord_map' (see below)
#
#-----------------------------------------------------------------------------------

class FETS2D4Qxfem(XFETSEval):
    
    parent_elem = FETS2D9Q()
        
    debug_on = True
    
    def ls_function(self, r,s): 
        return r - 0.1*s + 0.5#Y-0.2*X#temp
    # Dimensional mapping
    dim_slice = slice(0, 2)
    
    n_nodal_dofs = Int(4)
#    vtk_cell_type = 9 #Quad#
    vtk_cell_types = 'Triangle' #Triangle
    
    vtk_point_ip_map = [0,1,3,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
    
    def ls_fn_r(self,r,s):
        return self.ls_function(r,s)

    def ls_fn_s(self,s,r):
        return self.ls_function(r,s)

    pt_sets = Property(List,depends_on = 'ls_function')
    @cached_property
    def _get_pt_sets(self):
        inter_pts = self._get_inter_pts()
        pos_pts, neg_pts = self._get_ls_nodes()
        return [pos_pts, neg_pts, inter_pts]
        
    def _get_inter_pts(self):
        inter_pts = []
        for c_coord in [-1,1]:
            args = (c_coord)
            s_coord = self._get_intersect_pt(self.ls_fn_s, args)
            r_coord = self._get_intersect_pt(self.ls_fn_r, args)
            if s_coord:
                inter_pts.append([c_coord,s_coord])
            if r_coord:
                inter_pts.append([r_coord,c_coord])
        return array(inter_pts)
    
    def _get_intersect_pt(self, fn, args):
        try:
            return brentq(fn,-1,1,args=args)

        except ValueError:
            return
    
    def _get_ls_nodes(self):
        coords = array([[-1,-1],[1,-1],[1,1],[-1,1]])
        r,s = coords.T
        ls_fn = frompyfunc( self.ls_function, 2, 1 )
        ls_vals = ls_fn(r,s)
        pos_arr = coords[ls_vals[:] > 0.]
        neg_arr = coords[ls_vals[:] < 0.]
        return pos_arr, neg_arr
    
    
#    pt_sets =        [array([[-1., -1.],#temporary hack, later replaced by interst algorithm
#                            [ -1., 1.]]),#positive, negative, intersection pts
#                     array([[ 1., 1],
#                            [ 1., -1]]),
#                      array([[ 0.2, -1.],
#                            [0.2, 1.]])]
    
    node_ls_values = Property(Array, depends_on = 'parent_elem')
    @cached_property
    def _get_node_ls_values(self):
        X, Y = array(self.dof_r).T
        #print "coords " ,array(self.dof_r).T
        vect_fn = frompyfunc( self.ls_function, 2, 1 )
        values = vect_fn( X, Y )
        return sign(values)
    

    #---------------------------------------------------------------------
    # Method required to represent the element geometry
    #---------------------------------------------------------------------
 
    #---------------------------------------------------------------------
    # Method delivering the shape functions for the field variables and their derivatives
    #---------------------------------------------------------------------
    
    def get_N_mtx(self,r_pnt):
        '''
        Returns the matrix of the shape functions used for the field approximation
        containing zero entries. The number of rows corresponds to the number of nodal
        dofs. The matrix is evaluated for the specified local coordinate r.
        '''
        p_N_mtx = self.get_parent_N_mtx(r_pnt)
        #print "N parent ",p_N_mtx
        N_e_list = [p_N_mtx[:,i*2:i*2+2]* \
                (sign(self.ls_function(r_pnt[0],r_pnt[1]))-self.node_ls_values[i])\
                 for i in range(0,self.parent_nnodes)]#TODO:this is just 2D
        N_e_mtx = hstack(N_e_list)
        N_enr_mtx = zeros((2,self.parent_nnodes*2*2))#TODO:this is just 2D
        for i in range(0,self.parent_nnodes ):
            N_enr_mtx[:,i*4:i*4+2] = p_N_mtx[:,i*2:i*2+2]
            N_enr_mtx[:,i*4+2:i*4+4] = N_e_mtx[:,i*2:i*2+2] 
        return N_enr_mtx

    def get_dNr_psi_mtx(self,r_pnt):
        '''
        Return the derivatives of the shape functions
        '''
        p_dNr_mtx = self.get_parent_dNr_mtx(r_pnt)

        dNr_mtx = array( [ p_dNr_mtx[:,i]* \
                           (sign(self.ls_function(r_pnt[0],r_pnt[1]))\
                           - self.node_ls_values[i]) for i in range(0,self.parent_nnodes) ])
        return dNr_mtx.T
    
#----------------------- example --------------------

def example_with_new_domain():
    from ibvpy.api import \
        TStepper as TS, RTraceGraph, RTraceDomainListField, TLoop, \
        TLine, BCDof, IBVPSolve as IS
    from ibvpy.mats.mats2D.mats2D_elastic.mats2D_elastic import MATS2DElastic
    from ibvpy.api import BCDofGroup
    fets_eval = FETS2D4Qxfem(mats_eval = MATS2DElastic(E= 1.,nu=0.)) 
    

    from ibvpy.mesh.fe_grid import FEGrid
    from ibvpy.mesh.fe_level_set_domain import FELevelSetDomain
    from mathkit.mfn.mfn_line.mfn_line import MFnLineArray

    # Discretization
    domain = FEGrid( coord_max = (1.,1.,0.), 
                           shape   = (1, 1),
                           fets_eval = fets_eval)
    
 
    mf1 = MFnLineArray( #xdata = arange(10),
                       ydata = array([0.,0.1,0.2,0.3,0.3,]) )
    mf2 = MFnLineArray( #xdata = arange(10),
                       ydata = array([0., 0., 0., 0.,0.,0.1,0.2,0.3 ]) )
                                         
    tstepper = TS( 
         sdomain = domain,
         bcond_list =  [ BCDofGroup( var='u', value = 0., dims = [0],
                                  get_dof_method = domain.get_left_dofs ),   
#                         BCDofGroup( var='u', value = 0., dims = [1],
#                                  get_dof_method = domain.get_top_dofs ), 
                        BCDofGroup( var='u', value = 0., dims = [1],
                                  get_dof_method = domain.get_bottom_left_dofs ), 
                           BCDofGroup( var='u', value = 1., dims = [0],
                                       time_function = mf1.get_value,
                                  get_dof_method = domain.get_right_dofs ),                                                   
                         BCDofGroup( var='u', value = 1., dims = [1],
                                  time_function = mf2.get_value,
                                  get_dof_method = domain.get_right_dofs ) ],
         rtrace_list =  [ 
#                     RTraceGraph(name = 'Fi,right over u_right (iteration)' ,
#                               var_y = 'F_int', idx_y = right_dof,
#                               var_x = 'U_k', idx_x = right_dof,
#                               record_on = 'update'),
#                         RTraceDomainField(name = 'Stress' ,
#                         var = 'sig_app', idx = 0,
#                         record_on = 'update'),
                    RTraceDomainListField(name = 'Displacement' ,
                                    var = 'u', idx = 0,
                                    record_on = 'update'),
                     RTraceDomainListField(name = 'Strain' ,
                                    var = 'eps', idx = 0,
                                    #position = 'int_pnts',
                                    record_on = 'update',
                                    warp = True),
#                    RTraceDomainField(name = 'N0' ,
#                                      var = 'N_mtx', idx = 0,
#                                      record_on = 'update')
                ]             
            )
    
    #import cProfile
    # Add the time-loop control
    #global tloop
    tloop = TLoop( tstepper = tstepper,
                   tline = TLine( min = 0.0,  step = 1., max = 10.0 ) )
               
#    cProfile.run('tloop.eval()', 'tloop_prof' )
#    
#    import pstats
#    p = pstats.Stats('tloop_prof')
#    p.strip_dirs()
#    print 'cumulative'
#    p.sort_stats('cumulative').print_stats(20)
#    print 'time'
#    p.sort_stats('time').print_stats(20)
    
    tloop.eval()
    # Put the whole thing into the simulation-framework to map the
    # individual pieces of definition into the user interface.
    #
    from ibvpy.plugins.ibvpy_app import IBVPyApp
    app = IBVPyApp( ibv_resource = tloop )
    app.main()    
    
if __name__ == '__main__':
    example_with_new_domain()
