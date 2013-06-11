
from enthought.traits.api import \
     Array, Bool, Callable, Enum, Float, HasTraits, Interface, implements, \
     Instance, Int, Trait, Str, Enum, Callable, List, TraitDict, Any, \
     on_trait_change, Tuple, WeakRef, Delegate, Property, cached_property

from enthought.traits.ui.api import \
     Item, View, HGroup, ListEditor, VGroup, Group

from enthought.traits.ui.menu import \
     NoButtons, OKButton, CancelButton, Action, CloseAction, Menu, \
     MenuBar, Separator

from math  import \
     pow, fabs

from numpy import \
     array, zeros, int_, float_, ix_, dot, linspace, hstack, vstack, arange

from scipy.linalg import \
     inv, det

import time

from ibvpy.fets.fets_eval import FETSEval, RTraceEvalElemFieldVar

#-----------------------------------------------------------------------------
# FEQ4
#-----------------------------------------------------------------------------

class FEQ4T(FETSEval):

    debug_on = True
    
    # Material parameters
    #
    k = Float( 1.0, label = 'conductivity' )

    # Dimensional mapping
    #
    dim_slice = slice(0,2)

    # System mapping parameters
    #
    n_e_dofs = Int(4)
    n_nodal_dofs = Int(1)

    # Order of node positions for the formulation of shape function
    # (isoparametric formulation)
    _node_coord_map = Array( 'float_', (4,2), 
                             [[-1.,-1.],
                              [ 1.,-1.],
                              [ 1., 1.],
                              [-1., 1.]] )
    # Integration parameters
    #
    ngp_r = 2
    ngp_s = 2

    # Field visualization attributes
    #
    field_entity_type = 'quad'

#    # Corner nodes are used for visualization 
#    vtk_r = [[-1.,-1.],
#                        [ 1.,-1.],
#                        [ 1., 1.],
#                        [-1., 1.]]
#    field_faces = [[0,1,2,3]]

    # 4 corner nodes, 4 edge nodes and 1 interior nodes 
    vtk_r = [[-1.,-1.],
                        [ 0.,-1.],
                        [ 1.,-1.],
                        [-1., 0.],
                        [ 0., 0.],
                        [ 1., 0.],
                        [-1., 1.],
                        [ 0., 1.],
                        [ 1., 1.]]
    field_faces = [[0,1,4,3],
                   [1,2,5,4],
                   [3,4,7,6],
                   [4,5,8,7]]
    
    #---------------------------------------------------------------------
    # Method required to represent the element geometry
    #---------------------------------------------------------------------
    def get_N_geo_mtx(self, r_pnt):
        '''
        Return the value of shape functions for the specified local coordinate r
        '''
        cx = self._node_coord_map
        Nr = array( [[ 1/4.*(1 + r_pnt[0]*cx[i,0])*(1 + r_pnt[1]*cx[i,1]) for i in range(0,4) ]] )
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
        return self.get_N_geo_mtx(r_pnt)
        
    def get_dNr_mtx(self,r_pnt):
        '''
        Return the derivatives of the shape functions
        '''
        return self.get_dNr_geo_mtx(r_pnt)
    
    def get_B_mtx( self, r_pnt, X_mtx ):
        J_mtx = self.get_J_mtx(r_pnt,X_mtx)
        dNr_mtx = self.get_dNr_mtx( r_pnt )
        dNx_mtx = dot( inv( J_mtx ), dNr_mtx  )
        return dNx_mtx

    def get_mtrl_corr_pred(self, sctx, eps_mtx):
        D_mtx = array([[ self.k,0 ],[0,self.k]])
        sig_mtx = dot( D_mtx, eps_mtx )
        return sig_mtx, D_mtx

#    # @todo: (alex): is this method used?
#    def get_X_mtx( self, X_mtx ):
#        X_mtx = self._get_X_mtx()
#        # make sure that the supplied global coordinate is within the
#        # element domain
#        return X_mtx

#----------------------- example --------------------

if __name__ == '__main__':
    from core.api import \
        TStepper as TS, MGridDomain, RTraceGraph, RTraceDomainField, TLoop, \
        TLine, BCDof, IBVPSolve as IS, DOTSEval

    fets_eval = FEQ4T()

    single_elem = False
    if single_elem:
        from ibvp_solve.melem_domain import MElemDomain, RTraceElemField
        quad_elem = MElemDomain( X = [[-5,-5],
                                      [10, 0],
                                      [10,10],
                                      [0, 10]])
        
        ts = TS( tse = fets_eval,
                 sdomain = quad_elem,
                 bcond_list = [ BCDof(var='u', dof = i, value = 0.) for i in [0,3] ] +  
                           [ BCDof(var='f', dof = i, value = 10 ) for i in [1,2] ],
                 rtrace_list = [ RTraceElemField(name = 'Temperature' ,
                                          var = 'u', idx = 0,
                                          update_on = 'update'),
                             RTraceElemField(name = 'Flux' ,
                                          var = 'eps', idx = 0,
                                          update_on = 'update' ),
                             RTraceElemField(name = 'Jacobi determinant' ,
                                          var = 'J_det', idx = 0,
                                          update_on = 'update'),
                             RTraceElemField(name = 'Shape functions' ,
                                          var = 'N_mtx', idx = 0,
                                          update_on = 'update')
                                          ]
                 )
    
        # Add the time-loop control
        #
        tl = TLoop( tstepper = ts,
                 DT = 0.5,
                 tline  = TLine( min = 0.0,  max = 1.0 ))
    
        tl.eval()
        # Put the whole stuff into the simulation-framework to map the
        # individual pieces of definition into the user interface.
        #
        sim = IS( tloop = tl )
        sim.configure_traits()
 
    else:       
        # Tseval for a discretized line domain
        #
        tseval  = DOTSEval( fets_eval = FEQ4T() )
        
        # Discretization
        #
        line_domain = MGridDomain( lengths = (3.,1.,0.), 
                                      shape = (8,8,0), 
                                      n_e_nodes_geo = (1,1,0), 
                                      n_e_nodes_dof = (1,1,0), 
                                      node_map_geo = [0,1,3,2], 
                                      node_map_dof = [0,1,3,2] ) 
        
        # Put the tseval (time-stepper) into the spatial context of the
        # discretization and specify the response tracers to evaluate there.
        #
        right_dof = 2
        ts = TS( tse = tseval,
             sdomain = line_domain,
             bcond_list =  [ BCDof(var='u', dof = i, value = 0.) for i in line_domain.get_left_dofs()[:,0]] +  
                        [ BCDof(var='u', dof = i, value = 20 ) for i in line_domain.get_right_dofs()[:,0]],
             rtrace_list = [ RTraceGraph(name = 'Fi,right over u_right (iteration)' ,
                      var_y = 'F_int', idx_y = right_dof,
                      var_x = 'U_k', idx_x = right_dof,
                      update_on = 'update'),
#                      RTraceDomainField(name = 'Flux field' ,
#                      var = 'eps', idx = 0,
#                      update_on = 'update'),
                      RTraceDomainField(name = 'Temperature' ,
                      var = 'u', idx = 0,
                      update_on = 'update'),
#                             RTraceDomainField(name = 'Shape functions' ,
#                                          var = 'N_mtx', idx = 0,
#                                          update_on = 'update')
                      
                 ]             
             )
    
        # Add the time-loop control
        #
        tl = TLoop( tstepper = ts,
             DT = 0.5,
             tline  = TLine( min = 0.0,  max = 1.0 ))
    
        tl.eval()    
        # Put the whole stuff into the simulation-framework to map the
        # individual pieces of definition into the user interface.
        #
        sim = IS( tloop = tl )
        sim.configure_traits()