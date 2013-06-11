
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

class FETS2D4Q_disc(FETSEval):
    
    debug_on = True

    # Dimensional mapping
    dim_slice = slice(0, 2)
    
    n_e_dofs = Int(8)
    t = Float( 1.0, label = 'thickness' )
    E = Float( 1.0, label = "Young's modulus" )
    nu = Float( 0.2, label = "Poison's ratio" )

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
    vtk_point_ip_map = [0, 1, 3, 2]
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

if __name__ == '__main__':
    from ibvpy.api import \
        TStepper as TS, MGridDomain, RTraceGraph, RTraceDomainField, TLoop, \
        TLine, BCDof, IBVPSolve as IS, DOTSEval
    #from lib.mats.mats2D.mats_cmdm2D.mats_mdm2d import MACMDM
#    from lib.mats.mats2D.mats2D_sdamage.mats2D_sdamage import MATS2DScalarDamage
#    from lib.mats.mats2D.mats2D_sdamage.strain_norm2d import *
    from ibvpy.mats.mats2D.mats2D_elastic.mats2D_elastic import MATS2DElastic
#    fets_eval = FETS2D4Q12U(mats_eval = MATS2DScalarDamage(strain_norm = Euclidean())) 
    fets_eval = FETS2D4Q(mats_eval = MATS2DElastic()) 
    #fets_eval = FETS2D9Q(mats_eval = MACMDM())            
    # Tseval for a discretized line domain
    tseval  = DOTSEval( fets_eval = fets_eval )

    from ibvpy.mesh.mgrid_domain import MeshGridAdaptor

    # Define a mesh domain adaptor as a cached property to 
    # be constracted on demand
    mgrid_adaptor = MeshGridAdaptor( n_nodal_dofs = 2,
                                     # NOTE: the following properties must be defined and 
                                     # must correspond to the used element formulation
                                     n_e_nodes_geo = (1,1,0), 
                                     n_e_nodes_dof = (1,1,0), 
                                     node_map_geo = [0,1,3,2], 
                                     node_map_dof = [0,1,3,2] )

    # Discretization
    domain = MGridDomain( lengths = (3.,1.,0.), 
                             shape = (8,8,0), 
                             adaptor = mgrid_adaptor )
                                         
    right_dof = 2
    tstepper = TS( tse = tseval,
         sdomain = domain,
         # conversion to list (square brackets) is only necessary for slicing of 
         # single dofs, e.g "get_left_dofs()[0,1]"
         bcond_list =  [ BCDof(var='u', dof = i, value = 0.) for i in  domain.get_left_dofs()[:,0]  ] +
                    [ BCDof(var='u', dof = i, value = 0.) for i in [domain.get_left_dofs()[0,1]] ] +    
                    [ BCDof(var='u', dof = i, value = 0.002 ) for i in domain.get_right_dofs()[:,0] ],
         rtrace_list =  [ 
                     RTraceGraph(name = 'Fi,right over u_right (iteration)' ,
                               var_y = 'F_int', idx_y = right_dof,
                               var_x = 'U_k', idx_x = right_dof,
                               record_on = 'update'),
                         RTraceDomainField(name = 'Stress' ,
                         var = 'sig_app', idx = 0,
                         record_on = 'update'),
                     RTraceDomainField(name = 'Displacement' ,
                                    var = 'u', idx = 0,
                                    record_on = 'update',
                                    warp = True),
#                             RTraceDomainField(name = 'N0' ,
#                                          var = 'N_mtx', idx = 0,
#                                          record_on = 'update')
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
    