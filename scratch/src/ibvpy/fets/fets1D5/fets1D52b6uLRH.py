from enthought.traits.api import \
    Int, implements, List, Array, Property, cached_property

from ibvpy.fets.fets_eval import IFETSEval, FETSEval


from numpy import array, dot, identity, hstack

from scipy.linalg import \
     inv

#-----------------------------------------------------------------------------
# FEBar1D
#-----------------------------------------------------------------------------

class FETS1D52B6ULRH(FETSEval):
    '''
    Fe Bar 2 nodes, deformation
    '''
    
    implements( IFETSEval )
    
    debug_on = True
    
    # Dimensional mapping
    dim_slice = slice(0, 2)    
    
    n_e_dofs = Int(6)
    n_nodal_dofs = Int(1)

    dof_r =  [[-1,-1],[1,-1],[1,1],[-1,1],[0,-1],[0,1]]          
    geo_r =  [[-1,-1],[1,-1],[1,1],[-1,1]]
    
    #vtk_point_ip_map = [0, 2, 2, 0, 1, 1]
    #field_lines = [[0, 1]]

    vtk_r = [[-1.,-1.], [1.,-1.], [1.,1.], [-1.,1.],[0.,-1.], [0.,1.]]
    vtk_cells = [[0, 4, 1],[2, 5, 3],[0, 3],[4, 5],[1, 2]]
    vtk_cell_types = ['QuadraticEdge','QuadraticEdge','Line','Line','Line']
    
    def setup(self, sctx):
        '''
        overloading default setup - sending number of slip component 
        as additional parameter to materialmodel
        @param sctx:
        '''
        i = 0
        for gp in self.ip_coords:
            sctx.mats_state_array = sctx.elem_state_array[(i * self.m_arr_size): ((i+1)*self.m_arr_size)]
            sctx.slip_comp = 3
            self.mats_eval.setup( sctx )
            i += 1
            
    def get_mp_state_array_size( self, sctx ):
        '''Get the size of the state array for a single material point.
        overloading default setup - sending number of damage params
        as additional parameter to materialmodel
        '''
        return 3#hack for this combination

    def get_state_array_size( self ):
        return self.m_arr_size * self.ip_coords.shape[0] #TODO:synchronize with fets eval - this is more robust


    # Integration parameters
    #
    ip_coords = Property()#TODO: standard overloading does not work, has to be property?
    @cached_property
    def _get_ip_coords(self):
        return  array([[-1.,0., 0.],[0., 0., 0.],[1., 0., 0.]])
    
    ip_weights = Property()#TODO: standard overloading does not work, has to be property?
    @cached_property
    def _get_ip_weights(self):
        return array([[1./3.],[4./3.],[1./3.]], dtype = float)

    def get_N_geo_mtx( self, r_pnt ):
        '''
        Return geometric shape functions
        @param r_pnt:
        '''
        cx = array( self.geo_r, dtype = 'float_' )
        N_mtx = array( [[1/4.*(1 + r_pnt[0]*cx[i,0])*(1 + r_pnt[1]*cx[i,1]) 
                      for i in range(0,4) ]] )
        return N_mtx
    
 

    def get_N_mtx(self, r_pnt):
        '''
        Return shape functions
        @param r_pnt:local coordinates
        '''
        r = r_pnt[0]
        s = r_pnt[1]
        N1 = r/2.*(r+1.)
        N2 = r/2.*(r-1.)
        N3 = 1.-r*r
        Nf = (1.-s)/2.
        Nm = (1.+s)/2.
        N_mtx = array([[N2*Nf,N1*Nf,N1*Nm,N2*Nm,N3*Nf,N3*Nm]])
        return N_mtx

   
    def get_B_mtx( self, r_pnt, X ):
        '''
        Return kinematic matrix
        @param r:local coordinates
        @param X:global coordinates
        '''
        r = r_pnt[0]
        L = (X[1,0]-X[0,0])
        B1 = (2*r-1.)/L
        B2 = (2*r+1.)/L
        B3 = -4.*r/L
        B_mtx=array([[1., 0., 0., -1., 0., 0.],#s_l
                     [0., 1., -1., 0., 0., 0.],#s_r
                     [0., 0., 0., 0., 1.,-1.],#s_m
                     [B1, B2, 0., 0., B3, 0],#eps_f
                     [0., 0., B2, B1, 0., B3]])#eps_m
        return B_mtx
    
    def _get_J_det(self,r_pnt3d,X_mtx):
        return (X_mtx[1,0]- X_mtx[0,0])/2.

#----------------------- example --------------------

def example_with_new_domain():    
    from ibvpy.api import \
        TStepper as TS, RTraceGraph, RTraceDomainListField, TLoop, \
        TLine, IBVPSolve as IS, DOTSEval
    from ibvpy.api import BCDofGroup
    from ibvpy.mats.mats1D5.mats1D5bond_elastic_frictional import MATS1D5Bond
    from ibvpy.mesh.fe_grid import FEGrid
    from mathkit.mfn.mfn_line.mfn_line import MFnLineArray
        
    fets_eval = FETS1D52B6ULRH(mats_eval = MATS1D5Bond(Ef = 17000.,
                                                    Af = 2.65e-6/4.,
                                                    Am = 2.65e-6/4.,
                                                    Em = 17000.,
                                                    tau_max = 8.23 * 2,
                                                    tau_fr = 8.23  * 2 ,
                                                    s_cr = 0.030e-3 * 10 )) 
    # Discretization

    domain = FEGrid( coord_max = (1.,.1,0.), #new domain
                           shape   = (1,1),
                           fets_eval = fets_eval)
                                         
    ts = TS( dof_resultants = True,
         sdomain = domain,
         # conversion to list (square brackets) is only necessary for slicing of 
         # single dofs, e.g "get_left_dofs()[0,1]"
         bcond_list =  [BCDofGroup(var='u', value = 0.,dims = [0],
                               get_dof_method = domain.get_left_dofs),\
                      # imposed displacement for all right dofs in y-direction:
#                        BCDofGroup(var='u', value = 0., dims = [0],
#                            get_dof_method = domain.get_top_right_dofs ),
                        BCDofGroup(var='u', value = 1.e-3, dims = [0],
                            get_dof_method = domain.get_bottom_right_dofs )],
         rtrace_list =  [ RTraceGraph(name = 'Fi,right over u_right (iteration)' ,
                               var_y = 'F_int', idx_y = 1,
                               var_x = 'U_k', idx_x = 1),
                        RTraceDomainListField(name = 'Debonding' ,
                                var = 'debonding', idx = 0 ),
                        RTraceDomainListField(name = 'Displacement' ,
                                var = 'u', idx = 0),
#                             RTraceDomainListField(name = 'N0' ,
#                                          var = 'N_mtx', idx = 0,
#                                          record_on = 'update')
                      
                ]             
            )
    
    # Add the time-loop control
    tloop = TLoop( tstepper = ts,
             DT = 1.,
             tline  = TLine( min = 0.0,  max = 1.0 ))
    
    tloop.eval()    
    # Put the whole stuff into the simulation-framework to map the
    # individual pieces of definition into the user interface.
    #
    from ibvpy.plugins.ibvpy_app import IBVPyApp
    app = IBVPyApp( ibv_resource = tloop )
    app.main()

if __name__ == '__main__':
    example_with_new_domain()