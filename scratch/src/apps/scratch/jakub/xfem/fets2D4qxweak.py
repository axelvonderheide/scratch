
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
     identity, unique, average, frompyfunc, fabs, linalg, absolute

from scipy.linalg import \
     inv

from fets_eval_enr import FETSEvalEnr
from ibvpy.fets.fets2D.fets2D4q import FETS2D4Q

from fets_eval_xfem import FETSEvalXFEM

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

class FETS2D4Qxfem(FETSEvalXFEM):
    
    base_elem = FETS2D4Q()
        
    debug_on = True
    def ls_function(self, X,Y): 
        return X#Y-0.2*X#temp
    # Dimensional mapping
    dim_slice = slice(0, 2)
    
    n_e_dofs = Int(16)
    n_nodal_dofs = Int(4)
#    vtk_cell_type = 9 #Quad#
    vtk_cell_types = 'Triangle' #Triangle
    # Integration parameters
    #
    ngp_r = 2
    ngp_s = 2
    
    vtk_point_ip_map = [0,1,3,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]

    
    
    pt_sets =        [array([[-1., -1.],#temporary hack, later replaced by interst algorithm
                            [ -1., 1.]]),#positive, negative, intersection pts
                     array([[ 1., 1],
                            [ 1., -1]]),
                      array([[ 0., -1.],
                            [0., 1.]])]
    
    
    
#    vtk_r = [[-1., -1.],
#                        [ 1., -1.],
#                        [ 1., 1.],
#                        [-1., 1.],
#                        [-1., -0.01],
#                        [ 1., -0.01],
#                        [ -1., 0.01],
#                        [1., 0.01]] 
#   
#    field_faces = [[0, 1, 5, 4],[6,7,2,3]]
    

    #---------------------------------------------------------------------
    # Method required to represent the element geometry
    #---------------------------------------------------------------------
  
    #---------------------------------------------------------------------
    # Method delivering the shape functions for the field variables and their derivatives
    #---------------------------------------------------------------------
    def get_parent_N_mtx(self,r_pnt):
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
    
    def get_N_mtx(self,r_pnt):
        '''
        Returns the matrix of the shape functions used for the field approximation
        containing zero entries. The number of rows corresponds to the number of nodal
        dofs. The matrix is evaluated for the specified local coordinate r.
        '''
        p_N_mtx = self.get_parent_N_mtx(r_pnt)

        X, Y = array(self.dof_r).T
        vect_fn = frompyfunc( self.ls_function, 2, 1 )
        values = vect_fn( X, Y )
        sum_vals =  abs(sum(sum( p_N_mtx[:,i*2:i*2+2]* \
                           (absolute(self.ls_function(r_pnt[0],r_pnt[1]))\
                           - values[i]) for i in range(0,4) )))
        
        N_e_list = [p_N_mtx[:,i*2:i*2+2]* \
                (absolute(self.ls_function(r_pnt[0],r_pnt[1]))-sum_vals)\
                 for i in range(0,4)]
        N_e_mtx = hstack(N_e_list)
        N_enr_mtx = hstack((p_N_mtx[:,:2],N_e_mtx[:,:2],\
                            p_N_mtx[:,2:4],N_e_mtx[:,2:4],\
                            p_N_mtx[:,4:6],N_e_mtx[:,4:6],\
                            p_N_mtx[:,6:],N_e_mtx[:,6:]))#TODO: generalize that
        
        return N_enr_mtx

    def get_parent_dNr_mtx(self,r_pnt):
        '''
        Return the derivatives of the shape functions
        '''
        return self.get_dNr_geo_mtx(r_pnt)

    def get_dNr_psi_mtx(self,r_pnt):
        '''
        Return the derivatives of the shape functions
        '''
        p_dNr_mtx = self.get_parent_dNr_mtx(r_pnt)
        #print "parent dNr ",p_dNr_mtx
        X, Y = array(self.dof_r).T
        vect_fn = frompyfunc( self.ls_function, 2, 1 )
        values = vect_fn( X, Y )
        sum_vals =  abs(sum(sum(  p_dNr_mtx[:,i]* \
                           (absolute(self.ls_function(r_pnt[0],r_pnt[1]))\
                           - values[i]) for i in range(0,4) )))
       # print "sum vals ", sum_vals
        dNr_geo = array( [ p_dNr_mtx[:,i]* \
                           (absolute(self.ls_function(r_pnt[0],r_pnt[1]))\
                           - sum_vals) for i in range(0,4) ])
        return dNr_geo.T
    
    def get_B_mtx( self, r_pnt, X_mtx ):
        J_mtx = self.get_J_mtx(r_pnt,X_mtx)
        #print "loc_coords ", r_pnt
        dNr_mtx = self.get_parent_dNr_mtx( r_pnt )
        dNr_psi_mtx = self.get_dNr_psi_mtx( r_pnt )
        #print "dNr ",dNr_mtx
        #print "dNr_psi ",dNr_psi_mtx
        dNr_enr_mtx = hstack((dNr_mtx[:,:1],dNr_psi_mtx[:,:1],\
                              dNr_mtx[:,1:2],dNr_psi_mtx[:,1:2],\
                              dNr_mtx[:,2:3],dNr_psi_mtx[:,2:3],\
                              dNr_mtx[:,3:],dNr_psi_mtx[:,3:]))
        #print "dNr_both ",dNr_enr_mtx
        dNx_enr_mtx = dot( inv( J_mtx ), dNr_enr_mtx  )
        Bx_mtx = zeros( (3,16 ), dtype = 'float_' )
        for i in range(0,8):
            Bx_mtx[0,i*2]   = dNx_enr_mtx[0,i]
            Bx_mtx[1,i*2+1] = dNx_enr_mtx[1,i]
            Bx_mtx[2,i*2]   = dNx_enr_mtx[1,i]
            Bx_mtx[2,i*2+1] = dNx_enr_mtx[0,i]
        return Bx_mtx


    def get_mtrl_corr_pred(self, sctx, eps_mtx, d_eps, tn, tn1, eps_avg = None):
        r_pnt = sctx.r_pnt
        val = self.ls_function(r_pnt[0],r_pnt[1])
        if val >= 0.:
            sig_mtx, D_mtx = self.mats_eval.get_corr_pred(sctx, eps_mtx, d_eps, tn, tn1,)
        else:
            sig_mtx, D_mtx = self.mats_eval2.get_corr_pred(sctx, eps_mtx, d_eps, tn, tn1,) 
        return sig_mtx, D_mtx

#----------------------- example --------------------

def example_with_new_domain():
    from ibvpy.api import \
        TStepper as TS, RTraceGraph, RTraceDomainField, TLoop, \
        TLine, BCDof, IBVPSolve as IS, UniformDomainTSE
    from ibvpy.mats.mats2D.mats2D_elastic.mats2D_elastic import MATS2DElastic
    from ibvpy.api import BCDofGroup
    fets_eval = FETS2D4Qxfem(mats_eval = MATS2DElastic(E=5.,nu=0.2),
                             mats_eval2 = MATS2DElastic(E=2.,nu=0.2)) 
    tseval  = UniformDomainTSE( fets_eval = fets_eval,
                                cache_geo_matrices = True )

    from ibvpy.mesh.fe_grid import FEGrid
    from ibvpy.mesh.fe_level_set_domain import FELevelSetDomain
    from mathkit.mfn.mfn_line.mfn_line import MFnLineArray

    # Discretization
    domain = FEGrid( coord_max = (1.,1.,0.), 
                           shape   = (1, 1),
                           n_nodal_dofs = fets_eval.n_nodal_dofs,
                           dof_r = fets_eval.dof_r,
                           geo_r = fets_eval.geo_r)
    
    ls_domain = FELevelSetDomain( source_domain = domain,
                                  ls_function = lambda x,y: x**2 + y**2 - 1. )
 
    mf = MFnLineArray( #xdata = arange(10),
                       ydata = array([0,1,2,3]) )
                                         
    tstepper = TS( tse = tseval,
         sdomain = domain,
         bcond_list =  [ BCDofGroup( var='u', value = 0., dims = [0],
                                  get_dof_method = domain.get_left_dofs ),   
                         BCDofGroup( var='u', value = 0., dims = [1],
                                  get_dof_method = domain.get_bottom_left_dofs ), 
#                        BCDofGroup( var='u', value = 0., dims = [0],
#                                  get_dof_method = domain.get_top_left_dofs ),                                
                         BCDofGroup( var='u', value = 1., dims = [0],
                                  time_function = mf.get_value,
                                  get_dof_method = domain.get_right_dofs ) ],
         rtrace_list =  [ 
#                     RTraceGraph(name = 'Fi,right over u_right (iteration)' ,
#                               var_y = 'F_int', idx_y = right_dof,
#                               var_x = 'U_k', idx_x = right_dof,
#                               record_on = 'update'),
#                         RTraceDomainField(name = 'Stress' ,
#                         var = 'sig_app', idx = 0,
#                         record_on = 'update'),
                    RTraceDomainField(name = 'Displacement' ,
                                    var = 'u', idx = 0,
                                    record_on = 'update'),
                     RTraceDomainField(name = 'Strain' ,
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
                   tline = TLine( min = 0.0,  step = 0.2, max = 1.0 ) )
               
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
