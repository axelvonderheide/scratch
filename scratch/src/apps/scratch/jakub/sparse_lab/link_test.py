
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
     identity

from scipy.linalg import \
     inv

from ibvpy.fets.fets_eval import FETSEval

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

class FETS2D4Q(FETSEval):
    
    debug_on = True

    # Dimensional mapping
    dim_slice = slice(0, 2)
    
    n_e_dofs = Int(8)
    t = Float( 1.0, label = 'thickness' )
    
    vtk_cell_type = 9 #Quad
    # Integration parameters
    #
    ngp_r = 2
    ngp_s = 2

    field_entity_type = 'quad'
    # Corner nodes are used for visualization 
    vtk_r = [[-1., -1.],
                        [ 1., -1.],
                        [ 1., 1.],
                        [-1., 1.]] 
    #vtk_point_ip_map = [0, 1, 3, 2]
    field_faces = [[0, 1, 2, 3]]
    vtk_point_ip_map = [0,1,3,2]
    n_nodal_dofs = Int(2)
    
    # Order of node positions for the formulation of shape function
    #
    _node_coord_map = Array( 'float_', (4,2), 
                                 [[-1., -1.],
                                  [ 1., -1.],
                                  [ 1., 1.],
                                  [-1., 1.]] )

    #---------------------------------------------------------------------
    # Method required to represent the element geometry
    #---------------------------------------------------------------------
    def get_N_geo_mtx(self, r_pnt):
        '''
        Return the value of shape functions for the specified local coordinate r
        '''
        cx = self._node_coord_map
        Nr = array( [[1/4.*(1 + r_pnt[0]*cx[i,0])*(1 + r_pnt[1]*cx[i,1]) 
                      for i in range(0,4) ]] )
        return Nr

    def get_dNr_geo_mtx(self, r_pnt):
        '''
        Return the matrix of shape function derivatives.
        Used for the conrcution of the Jacobi matrix.

        @TODO - the B matrix is used
        just for uniaxial bar here with a trivial differential
        operator.
        '''
        cx = self._node_coord_map
        dNr_geo = array( [[ 1/4.*cx[i,0]*(1 + r_pnt[1]*cx[i,1]) for i in range(0,4) ],
                          [ 1/4.*cx[i,1]*(1 + r_pnt[0]*cx[i,0]) for i in range(0,4) ]])        
        return dNr_geo

    #---------------------------------------------------------------------
    # Method delivering the shape functions for the field variables and their derivatives
    #---------------------------------------------------------------------
    def get_N_mtx(self,r_pnt):
        '''
        Returns the matrix of the shape functions used for the field approximation
        containing zero entries. The number of rows corresponds to the number of nodal
        dofs. The matrix is evaluated for the specified local coordinate r.
        '''
        Nr_geo = self.get_N_geo_mtx(r_pnt)
        I_mtx = identity(self.n_nodal_dofs, float)
        N_mtx_list = [I_mtx*Nr_geo[0,i] for i in range(0,Nr_geo.shape[1])]
        N_mtx = hstack(N_mtx_list)
        return N_mtx

    def get_dNr_mtx(self,r_pnt):
        '''
        Return the derivatives of the shape functions
        '''
        return self.get_dNr_geo_mtx(r_pnt)
    
    def get_B_mtx( self, r_pnt, X_mtx ):
        J_mtx = self.get_J_mtx(r_pnt,X_mtx)
        dNr_mtx = self.get_dNr_mtx( r_pnt )
        dNx_mtx = dot( inv( J_mtx ), dNr_mtx  )
        Bx_mtx = zeros( (3,8 ), dtype = 'float_' )
        for i in range(0,4):
            Bx_mtx[0,i*2]   = dNx_mtx[0,i]
            Bx_mtx[1,i*2+1] = dNx_mtx[1,i]
            Bx_mtx[2,i*2]   = dNx_mtx[1,i]
            Bx_mtx[2,i*2+1] = dNx_mtx[0,i]
        return Bx_mtx


#----------------------- example --------------------

def example_with_new_domain():
    from ibvpy.api import \
        TStepper as TS, MGridDomain, RTraceGraph, RTraceDomainField, TLoop, \
        TLine, BCDof, IBVPSolve as IS, DOTSEval
    from ibvpy.mats.mats2D.mats2D_elastic.mats2D_elastic import MATS2DElastic
    from ibvpy.api import BCDofGroup
    fets_eval = FETS2D4Q(mats_eval = MATS2DElastic()) 
    tseval  = DOTSEval( fets_eval = fets_eval )

    from ibvpy.mesh.fe_grid import FEGrid
    from mathkit.mfn.mfn_line.mfn_line import MFnLineArray

    # Discretization
    domain = FEGrid( coord_max = (1.,1.,0.), 
                           shape   = (1, 2),
                           n_nodal_dofs = 2,
                           dof_r = [[-1,-1],[1,-1],[1,1],[-1,1]],
                           geo_r = [[-1,-1],[1,-1],[1,1],[-1,1]])
                         
#    bcg = BCDofGroup( var='u', value = 0., dims = [0],
#                   get_dof_method = domain.get_left_dofs )
#    bcg.setup( None )
#    print 'labels', bcg._get_labels()
#    print 'points', bcg._get_mvpoints()
#    
    mf = MFnLineArray( #xdata = arange(10),
                       ydata = array([0,1,2,3]) )
                                         
    right_dof = 2
    tstepper = TS( tse = tseval,
         sdomain = domain,
         bcond_list =  [ BCDofGroup( var='u', value = 0., dims = [0,1],
                                  get_dof_method = domain.get_bottom_dofs ), 
                        BCDofGroup( var='u', value = 0.2, dims = [0],
                                  # link_dofs = [6],
                                  #link_coeffs = [-1.],
                                  get_dof_method = domain.get_top_right_dofs ),          
                         BCDofGroup( var='f', value = 0.0, dims = [1],
                                  link_dofs = [5],
                                  link_coeffs = [1.],
                                  #time_function = mf.get_value,
                                  get_dof_method = domain.get_top_right_dofs ) , 
#                        BCDofGroup( var='f', value = 0.0, dims = [1],
#                                  link_dofs = [9],
#                                  link_coeffs = [1.],
#                                  #time_function = mf.get_value,
#                                  get_dof_method = domain.get_left_middle_dofs ) , 
                           ],
         rtrace_list =  [ 
#                     RTraceGraph(name = 'Fi,right over u_right (iteration)' ,
#                               var_y = 'F_int', idx_y = right_dof,
#                               var_x = 'U_k', idx_x = right_dof,
#                               record_on = 'update'),
#                         RTraceDomainField(name = 'Stress' ,
#                         var = 'sig_app', idx = 0,
#                         record_on = 'update'),
                     RTraceDomainField(name = 'Displacement' ,
                                    var = 'u', idx = 1,
                                    record_on = 'update',
                                    warp = True),
#                    RTraceDomainField(name = 'N0' ,
#                                      var = 'N_mtx', idx = 0,
#                                      record_on = 'update')
                ]             
            )
    
    # Add the time-loop control
    tloop = TLoop( tstepper = tstepper,
                   tline = TLine( min = 0.0,  step = 0.5, max = 1.0 ) )
               
    tloop.eval()

    # Put the whole thing into the simulation-framework to map the
    # individual pieces of definition into the user interface.
    #
    from ibvpy.plugins.ibvpy_app import IBVPyApp
    app = IBVPyApp( ibv_resource = tloop )
    app.main()    
    
if __name__ == '__main__':
    example_with_new_domain()
