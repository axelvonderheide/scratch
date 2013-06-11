'''
Created on Jun 14, 2009

@author: jakub
'''

from enthought.traits.api import \
     Array, Bool, Callable, Enum, Float, HasTraits, Interface, implements, \
     Instance, Int, Trait, Str, Enum, Callable, List, TraitDict, Any, Dict,\
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
     
from ibvpy.fets.fets_eval import RTraceEvalElemFieldVar
     
from ibvpy.fets.fets2D.fets2D4q import FETS2D4Q
from ibvpy.fets.fets2D.fets2D4q8u import FETS2D4Q8U
from ibvpy.fets.fets2D.fets2D4q9u import FETS2D4Q9U
from ibvpy.fets.fets2D.fets2D9q import FETS2D9Q

from ibvpy.fets.fets_ls.fets_ls_eval import FETSLSEval

class FETSCrackTF(FETSLSEval):

    x_slice = slice(0, 2)

    n_nodal_dofs = Property(Int, depends_on = 'parent_fets.n_nodal_dofs')
    @cached_property
    def _get_n_nodal_dofs(self):
        return self.parent_fets.n_nodal_dofs * 1.5
    
    n_e_dofs = Property(Int, depends_on = 'parent_fets.n_e_dofs')
    @cached_property
    def _get_n_e_dofs(self):
        return self.parent_fets.n_e_dofs * 1.5
 
               
               
    #---------------------------------------------------------------------
    # Method required to represent the element geometry
    #---------------------------------------------------------------------
 
    #---------------------------------------------------------------------
    # Method delivering the shape functions for the field variables and their derivatives
    #---------------------------------------------------------------------
    
    def get_N_mtx(self,r_pnt, node_ls_values, r_ls_value):
        '''
        Returns the matrix of the shape functions used for the field approximation
        containing zero entries. The number of rows corresponds to the number of nodal
        dofs. The matrix is evaluated for the specified local coordinate r.
        '''
        p_N_mtx = self.parent_fets.get_N_mtx(r_pnt)
        n_nodes = self.parent_fets.n_e_dofs/self.parent_fets.n_nodal_dofs
        p_nodal_dofs =  4
        p_value = sign(r_ls_value)
        N_e_mtx = hstack([p_N_mtx[:,i*4:i*4+2]* \
                    (p_value-sign(node_ls_values[i]))/2.\
                    for i in range(0,n_nodes)])

        N_enr_mtx = hstack((p_N_mtx,N_e_mtx))
# below version for the alternating standard and enriched dofs
#        N_enr_mtx = zeros((2,self.parent_fets.n_e_dofs *2))
#        for i in range(0,self.parent_fets.n_e_dofs ):
#            N_enr_mtx[:,i*4:i*4+2] = p_N_mtx[:,i*2:i*2+2]
#            N_enr_mtx[:,i*4+2:i*4+4] = N_e_mtx[:,i*2:i*2+2]
        return N_enr_mtx

    def get_dNr_psi_mtx(self,r_pnt, node_ls_values, r_ls_value):
        '''
        Return the derivatives of the shape functions
        '''
        p_dNr_mtx = self.get_dNr_mtx(r_pnt)
        p_value = sign(r_ls_value)
        n_nodes = self.parent_fets.n_e_dofs/self.parent_fets.n_nodal_dofs
        dNr_mtx = hstack( [ p_dNr_mtx[:,i:i+1]* \
                           (p_value\
                           - sign(node_ls_values[i]))/2. for i in range(0,n_nodes) ])
        return dNr_mtx
    
    def get_B_mtx( self, r_pnt, X_mtx, node_ls_values, r_ls_value ):
        J_mtx = self.get_J_mtx(r_pnt,X_mtx)
        dNr_mtx = self.get_dNr_mtx( r_pnt )
        dNr_psi_mtx = self.get_dNr_psi_mtx( r_pnt, node_ls_values, r_ls_value )
        dNx_mtx = dot( inv( J_mtx ), dNr_mtx  )
        dNx_psi_mtx = dot( inv( J_mtx ), dNr_psi_mtx  )
        N_mtx = self.get_N_mtx(r_pnt, node_ls_values, r_ls_value)
        N_mtx_red = N_mtx[:2]-N_mtx[2:]
        n_nodes = self.parent_fets.n_e_dofs/self.parent_fets.n_nodal_dofs
        if dNr_mtx.shape[0] == 1:#1D#TODO:
            Bx_mtx = zeros( (1,  self.n_e_dofs), dtype = 'float_' )
            for i in range(0,n_nodes):
                Bx_mtx[0,i]   = dNx_mtx[0,i]               
                Bx_mtx[0,i+2] = dNx_psi_mtx[0,i]
 
        elif dNr_mtx.shape[0] == 2:#2D
            Bx_mtx = zeros( (8,  self.n_e_dofs), dtype = 'float_' )
            for i in range(0,n_nodes):
                
                Bx_mtx[0,i*4]   = dNx_mtx[0,i]
                Bx_mtx[1,i*4+1] = dNx_mtx[1,i]
                Bx_mtx[2,i*4]   = dNx_mtx[1,i]
                Bx_mtx[2,i*4+1] = dNx_mtx[0,i]
            
                Bx_mtx[3,i*4+2] = dNx_mtx[0,i]
                Bx_mtx[4,i*4+3] = dNx_mtx[1,i]
                Bx_mtx[5,i*4+2] = dNx_mtx[1,i]
                Bx_mtx[5,i*4+3] = dNx_mtx[0,i]
               
                Bx_mtx[0,i*2+self.parent_fets.n_e_dofs] = dNx_psi_mtx[0,i]
                Bx_mtx[1,i*2+self.parent_fets.n_e_dofs+1] = dNx_psi_mtx[1,i]
                Bx_mtx[2,i*2+self.parent_fets.n_e_dofs] = dNx_psi_mtx[1,i]
                Bx_mtx[2,i*2+self.parent_fets.n_e_dofs+1] = dNx_psi_mtx[0,i] 
                
            Bx_mtx[6:] = N_mtx_red#TODO:this is just 2d
        return Bx_mtx
    
    def map_eps(self, sctx, u):
        e_id = sctx.e_id
        p_id = sctx.p_id
        X_mtx = sctx.X
        r_pnt = sctx.loc
        B_mtx = self.get_B_mtx(r_pnt, X_mtx,
                               sctx.dots.dof_node_ls_values[e_id],
                               sctx.dots.vtk_ls_values[e_id][p_id] )
        return dot( B_mtx, u )
    
    
    def get_eps_m(self, sctx, u):
        e_id = sctx.e_id
        p_id = sctx.p_id
        X_mtx = sctx.X
        r_pnt = sctx.loc
        B_mtx = self.get_B_mtx(r_pnt, X_mtx,
                               sctx.dots.dof_node_ls_values[e_id],
                               sctx.dots.vtk_ls_values[e_id][p_id] )
        eps = dot( B_mtx, u )
        return array([[eps[0],eps[2]],[eps[2],eps[1]]])
    
    def get_eps_f(self, sctx, u):
        e_id = sctx.e_id
        p_id = sctx.p_id
        X_mtx = sctx.X
        r_pnt = sctx.loc
        B_mtx = self.get_B_mtx(r_pnt, X_mtx,
                               sctx.dots.dof_node_ls_values[e_id],
                               sctx.dots.vtk_ls_values[e_id][p_id] )
        eps = dot( B_mtx, u )
        return array([[eps[3],eps[5]],[eps[5],eps[4]]])
    
    def get_eps_b(self, sctx, u):
        e_id = sctx.e_id
        p_id = sctx.p_id
        X_mtx = sctx.X
        r_pnt = sctx.loc
        B_mtx = self.get_B_mtx(r_pnt, X_mtx,
                               sctx.dots.dof_node_ls_values[e_id],
                               sctx.dots.vtk_ls_values[e_id][p_id] )
        eps = dot( B_mtx, u )
        return array([[eps[6],0.],[0.,eps[7]]])
    
    def get_u_m(self, sctx, u):
        e_id = sctx.e_id
        p_id = sctx.p_id
        N_mtx = self.get_N_mtx( sctx.loc,
                                sctx.dots.dof_node_ls_values[e_id],
                                sctx.dots.vtk_ls_values[e_id][p_id] )
        return dot( N_mtx, u )[:2]
    
    def get_u_f(self, sctx, u):
        e_id = sctx.e_id
        p_id = sctx.p_id
        N_mtx = self.get_N_mtx( sctx.loc,
                                sctx.dots.dof_node_ls_values[e_id],
                                sctx.dots.vtk_ls_values[e_id][p_id] )
        return dot( N_mtx, u )[2:]
    
    rte_dict = Trait( Dict )
    def _rte_dict_default(self):
        '''
        RTraceEval dictionary with standard field variables.
        '''
        rte_dict = self._debug_rte_dict()
        for key, v_eval in self.mats_eval.rte_dict.items():

            # add the eval into the loop.
            #
            rte_dict[ key ] = RTraceEvalElemFieldVar( name = key, 
                                                      u_mapping = self.map_eps,
                                                      eval = v_eval )

        rte_dict.update( {'eps_m' : RTraceEvalElemFieldVar( eval = self.get_eps_m ), 
                          'u_m'   : RTraceEvalElemFieldVar( eval = self.get_u_m ),
                          'eps_f' : RTraceEvalElemFieldVar( eval = self.get_eps_f ), 
                          'u_f'   : RTraceEvalElemFieldVar( eval = self.get_u_f ),
                          'eps_b' : RTraceEvalElemFieldVar( eval = self.get_eps_b )} )

        return rte_dict
