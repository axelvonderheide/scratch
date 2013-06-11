
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
     identity, repeat, diag

from scipy.linalg import \
     inv

from ibvpy.fets.fets_eval import FETSEval, RTraceEvalElemFieldVar
from ibvpy.fets.fets2D.fets2D4q import FETS2D4Q

#-----------------------------------------------------------------------------------
# FETS2D4Q - 4 nodes iso-parametric quadrilateral element (2D, linear, Lagrange family)    
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

class FETS2D4Qtf(FETS2D4Q):
    
    debug_on = True

    # Dimensional mapping
    dim_slice = slice(0, 2)
    
    n_e_dofs = Int(24)
    
    # Integration parameters
    #
    ngp_r = 2
    ngp_s = 2

    # Corner nodes are used for visualization 
    vtk_r = [[-1., -1.],
                  [ 1., -1.],
                  [ 1.,  1.],
                  [-1.,  1.]] 
    vtk_cells = [[0, 1, 2, 3]]
    vtk_cell_types = 'Quad'
    
    #vtk_point_ip_map = [0, 1, 3, 2]
    vtk_point_ip_map = [0,1,3,2]
    n_nodal_dofs = Int(6)
    
    # Order of node positions for the formulation of shape function
    #
    dof_r = [[-1,-1],[1,-1],[1,1],[-1,1]]  
    geo_r = [[-1,-1],[1,-1],[1,1],[-1,1]]

    #---------------------------------------------------------------------
    # Method delivering the shape functions for the field variables and their derivatives
    #---------------------------------------------------------------------
    def get_N_s_mtx(self,r_pnt):
        '''
        Returns the matrix of the shape functions used for the field approximation
        containing zero entries. The number of rows corresponds to the number of nodal
        dofs. The matrix is evaluated for the specified local coordinate r.
        '''
        Nr_geo = self.get_N_geo_mtx(r_pnt)
        I_mtx = identity(2, float)
        N_mtx_list = [I_mtx*Nr_geo[0,i] for i in range(0,Nr_geo.shape[1])]
        N_mtx = hstack(N_mtx_list)
        return N_mtx    
    #---------------------------------------------------------------------
    # Method delivering the shape functions for the field variables and their derivatives
    #---------------------------------------------------------------------
    def get_N_mtx(self,r_pnt):
        '''
        Returns the matrix of the shape functions used for the field approximation
        containing zero entries. The number of rows corresponds to the number of nodal
        dofs. The matrix is evaluated for the specified local coordinate r.
        '''
        N_s_mtx = self.get_N_s_mtx(r_pnt)
        
        Nr_geo = self.get_N_geo_mtx(r_pnt)
        I_mtx = identity(4, float)
        N_mtx_list = [I_mtx*Nr_geo[0,i] for i in range(0,Nr_geo.shape[1])]
        N_e_mtx = hstack(N_mtx_list)
        N_mtx = zeros( (4, 24), dtype = 'float_' )
        for i in range(0,4):
            N_mtx[:,i*6:i*6+6] = hstack((vstack((N_s_mtx[:,i*2:i*2+2],
                                                 N_s_mtx[:,i*2:i*2+2])),
                                                 N_e_mtx[:,i*4:i*4+4]))
        return N_mtx                             


    def get_B_s_mtx( self, r_pnt, X_mtx ):
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
    
    def get_B_mtx( self, r_pnt, X_mtx ):
        B_s_mtx = self.get_B_s_mtx(r_pnt, X_mtx)
        
        J_mtx = self.get_J_mtx(r_pnt,X_mtx)
        dNr_mtx = self.get_dNr_mtx( r_pnt )
        dNx_mtx = dot( inv( J_mtx ), dNr_mtx  )
        Bx_mtx = zeros( (11,24 ), dtype = 'float_' )#9+2 strain component to 24 DOFs
        for i in range(0,4):
            Bx_mtx[0:3,i*6:i*6+2] = B_s_mtx[:,i*2:i*2+2]
            
            Bx_mtx[3,i*6+2] = dNx_mtx[0,i]
            Bx_mtx[4,i*6+3] = dNx_mtx[1,i]
            Bx_mtx[5,i*6+2] = dNx_mtx[1,i]
            Bx_mtx[5,i*6+3] = dNx_mtx[0,i]
            Bx_mtx[6,i*6+4] = dNx_mtx[0,i]
            Bx_mtx[7,i*6+5] = dNx_mtx[1,i]
            Bx_mtx[8,i*6+4] = dNx_mtx[1,i]
            Bx_mtx[8,i*6+5] = dNx_mtx[0,i]
            Bx_mtx[9,i*6:(i+1)*6] = [0.,0.,1.,0.,-1.,0.]
            Bx_mtx[10,i*6:(i+1)*6]= [0.,0.,0.,1.,0.,-1.]
        
        #print "B_mtx ", Bx_mtx
        return Bx_mtx

    def get_u_m(self, sctx, u):
        N_mtx = self.get_N_mtx( sctx.loc )
        return dot( N_mtx, u )[0:2]

    def get_u_f(self, sctx, u):
        N_mtx = self.get_N_mtx( sctx.loc )
        return dot( N_mtx, u )[2:4]

    def _rte_dict_default(self):
        '''
        RTraceEval dictionary with standard field variables.
        '''
        rte_dict = super( FETS2D4Qtf, self )._rte_dict_default()
        rte_dict['u_m'] = RTraceEvalElemFieldVar( eval = self.get_u_m,   ts = self )
        rte_dict['u_f'] = RTraceEvalElemFieldVar( eval = self.get_u_f,   ts = self )
        return rte_dict

#----------------------- example --------------------

def example_with_new_domain():
    from ibvpy.api import \
        TStepper as TS, RTraceGraph, RTraceDomainField, TLoop, \
        TLine, BCDof, IBVPSolve as IS, DOTSEval
    from ibvpy.mats.mats2D.mats2D_elastic.mats2D_elastic import MATS2DElastic
    from ibvpy.mats.mats2D.mats2D_sdamage.mats2D_sdamage import MATS2DScalarDamage
    from mats2D5bond import MATS2D5Bond

    from ibvpy.api import BCDofGroup
    fets_eval = FETS2D4Qtf(mats_eval = MATS2D5Bond()) 
    #fets_eval = FETS2D4Qtf(mats_eval = MATS2DScalarDamage()) 

    #print fets_eval.vtk_cell_data
    
    tseval  = DOTSEval( fets_eval = fets_eval,
                                cache_geo_matrices = True )

    from ibvpy.mesh.fe_grid import FEGrid
    from mathkit.mfn.mfn_line.mfn_line import MFnLineArray

    # Discretization
    domain = FEGrid( coord_max = (1.,2.,0.), 
                           shape   = (5, 10),
                           n_nodal_dofs = fets_eval.n_nodal_dofs,
                           dof_r = fets_eval.dof_r,
                           geo_r = fets_eval.geo_r )
                         
    mf = MFnLineArray( #xdata = arange(10),
                       ydata = array([0,1,2,3]) )
                                         
    tstepper = TS( tse = tseval,
         sdomain = domain,
         bcond_list =  [ BCDofGroup( var='u', value = 0., dims = [0,1,2,3],
                                  get_dof_method = domain.get_bottom_dofs ),
                         #BCDofGroup( var='u', value = 0., dims = [1],
                         #         get_dof_method = domain.get_bottom_dofs ),                                  
                         BCDofGroup( var='u', value = 1., dims = [4],
                                  #time_function = mf.get_value,
                                  get_dof_method = domain.get_top_right_dofs ) ],
         rtrace_list =  [ 
#                     RTraceGraph(name = 'Fi,right over u_right (iteration)' ,
#                               var_y = 'F_int', idx_y = right_dof,
#                               var_x = 'U_k', idx_x = right_dof,
#                               record_on = 'update'),
                        RTraceDomainField(name = 'Stress' ,
                        var = 'sig_app', idx = 0,
                        record_on = 'update'),
                        RTraceDomainField(name = 'Shear' ,
                        var = 'shear', idx = 0,
                        record_on = 'update'),
                     RTraceDomainField(name = 'Displ f' ,
                                    var = 'u_f', idx = 0,
                                    record_on = 'update',
                                    warp = True),
                    RTraceDomainField(name = 'Displ m' ,
                                    var = 'u_m', idx = 0,
                                    record_on = 'update',
                                    warp = True),
                    RTraceDomainField(name = 'Displ s' ,
                                    var = 'u', idx = 0,                                    record_on = 'update',
                                    warp = True),                               
#                    RTraceDomainField(name = 'N0' ,
#                                      var = 'N_mtx', idx = 0,
#                                      record_on = 'update')
                ]             
            )
    
    # Add the time-loop control
    #global tloop
    tloop = TLoop( tstepper = tstepper, #debug = True,
                   tline = TLine( min = 0.0,  step = .5, max = 1.0 ) )

    #import cProfile
    #cProfile.run('tloop.eval()', 'tloop_prof' )
    
    tloop.eval()
    
#    import pstats
#    p = pstats.Stats('tloop_prof')
#    p.strip_dirs()
#    print 'cumulative'
#    p.sort_stats('cumulative').print_stats(20)
#    print 'time'
#    p.sort_stats('time').print_stats(20)

    # Put the whole thing into the simulation-framework to map the
    # individual pieces of definition into the user interface.
    #
    from ibvpy.plugins.ibvpy_app import IBVPyApp
    app = IBVPyApp( ibv_resource = tloop )
    app.main()    
    
if __name__ == '__main__':
    example_with_new_domain()
