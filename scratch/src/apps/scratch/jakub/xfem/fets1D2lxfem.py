from enthought.traits.api import \
    Int, implements

from ibvpy.fets.fets_eval import IFETSEval, FETSEval


from numpy import array, dot, frompyfunc, hstack

from scipy.linalg import \
     inv
from ibvpy.fets.fets1D.fets1D2l import FETS1D2L

def signum(a):
    if a == 0.:
        return 0
    else:
        return a / abs(a)

#-----------------------------------------------------------------------------
# FEBar1D
#-----------------------------------------------------------------------------

class FETS1D2Lxfem(FETS1D2L):
    '''
    Fe Bar 2 nodes, deformation
    '''
    
    implements( IFETSEval )
    
    debug_on = True
    
    def ls_function(self, X): 
        return X-0.#temp
    # Dimensional mapping
    dim_slice = slice(0, 1)    
    
    n_e_dofs = Int(4)
    n_nodal_dofs = Int(2)

    dof_r = [[-1],[1]]
    geo_r = [[-1],[1]]

    vtk_r = [[-1.],[-0.1],[0.1],[1.]]
    vtk_cells = [[0, 1],[2,3]]


    vtk_point_ip_map = [0, 0, 0, 0]
    
    # Integration parameters
    #
    ngp_r = 2

    def get_N_mtx(self, r_pnt):
        '''
        Return shape functions
        @param r_pnt:local coordinates
        '''
        X = array([[-1],[1]])
        vect_fn = frompyfunc( self.ls_function, 1, 1 )
        values = vect_fn( X )
        sig_vals = []
        for val in values:
            sig_vals.append(signum(val))
        
        N_mtx = self.get_N_geo_mtx(r_pnt)
        
        N_e_mtx_list = [N_mtx[0,i]* \
                (signum(self.ls_function(r_pnt[0]))-sig_vals[i])\
                 for i in range(0,2)]
        N_e_mtx = array(N_e_mtx_list).T
        
        #print "N_mtx ",N_mtx
        #print "N_e_mtx ", N_e_mtx
        
        N_enr_mtx = hstack((N_mtx[:,:1],N_e_mtx[:,:1],\
                            N_mtx[:,1:2],N_e_mtx[:,1:2]))
        #print "N_enr_mtx ", N_enr_mtx
        return N_enr_mtx
        
        return self.get_N_geo_mtx(r_pnt)

    def get_dNr_psi_mtx(self, r_pnt):
        '''
        Return the derivatives of the shape functions
        '''
        return self.get_dNr_geo_mtx(r_pnt)
    
    def get_B_mtx( self, r, X ):
        '''
        Return kinematic matrix
        @param r:local coordinates
        @param X:global coordinates
        '''
        J_mtx = self.get_J_mtx(r, X)
        dNr_mtx = self.get_dNr_mtx( r )
        dNr_psi_mtx = self.get_dNr_psi_mtx(r)
        #print "dNr ",dNr_mtx
        #print "dNr_psi ",dNr_psi_mtx
        dNr_enr_mtx = hstack((dNr_mtx[:,:1],dNr_psi_mtx[:,:1],\
                              dNr_mtx[:,1:2],dNr_psi_mtx[:,1:2]))
        #print "dNr_both ",dNr_enr_mtx
        B_mtx = dot( inv( J_mtx ), dNr_enr_mtx  )
        return B_mtx

#----------------------- example --------------------

def example_with_new_domain():    
    from ibvpy.api import \
        TStepper as TS, RTraceGraph, RTraceDomainField, TLoop, \
        TLine, BCDofGroup, IBVPSolve as IS, DOTSEval
    from ibvpy.mats.mats1D.mats1D_elastic.mats1D_elastic import MATS1DElastic
    
    fets_eval = FETS1D2Lxfem(mats_eval = MATS1DElastic(E=10., A=1.))        

    # Tseval for a discretized line domain
    tseval  = DOTSEval( fets_eval = fets_eval )

    from ibvpy.mesh.fe_grid import FEGrid

    # Discretization
    domain = FEGrid( coord_max = (1.,0.,0.), 
                           shape   = (1,),
                           n_nodal_dofs = fets_eval.n_nodal_dofs,
                           dof_r = fets_eval.dof_r,
                           geo_r = fets_eval.geo_r)
                                                 
    ts = TS( tse = tseval,
             dof_resultants = True,
             sdomain = domain,
            bcond_list =  [ BCDofGroup( var='u', value = 0., dims = [0],
                                  get_dof_method = domain.get_left_dofs ),                                
                         BCDofGroup( var='u', value = 1., dims = [0],
                                  get_dof_method = domain.get_right_dofs ) ],
         rtrace_list =  [ RTraceGraph(name = 'Fi,right over u_right (iteration)' ,
                               var_y = 'F_int', idx_y = 0,
                               var_x = 'U_k', idx_x = 1),
                    RTraceDomainField(name = 'Stress' ,
                         var = 'sig_app', idx = 0),
                     RTraceDomainField(name = 'Displacement' ,
                                    var = 'u', idx = 0)
                      
                ]             
            )
    
    # Add the time-loop control
    tloop = TLoop( tstepper = ts,
                   tline  = TLine( min = 0.0,  step = 1, max = 1.0 ))
    
    print '---- result ----'
    print tloop.eval()
    print ts.F_int
    print ts.rtrace_list[0].trace.ydata
    
    # Put the whole stuff into the simulation-framework to map the
    # individual pieces of definition into the user interface.
    #
    from ibvpy.plugins.ibvpy_app import IBVPyApp
    app = IBVPyApp( ibv_resource = tloop )
    app.main()


if __name__ == '__main__':
    example_with_new_domain()