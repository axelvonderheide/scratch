'''
Created on Jun 14, 2009

@author: jakub
'''

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
     
from ibvpy.fets.fets2D.fets2D4q import FETS2D4Q
from ibvpy.fets.fets2D.fets2D4q8u import FETS2D4Q8U
from ibvpy.fets.fets2D.fets2D4q9u import FETS2D4Q9U
from ibvpy.fets.fets2D.fets2D9q import FETS2D9Q

from ibvpy.fets.fets_ls.fets_ls_eval import FETSLSEval

class FETSCrackRough(FETSLSEval):
    
    x_slice = slice(0,1)

    n_nodal_dofs = Property(Int, depends_on = 'parent_fets.n_nodal_dofs')
    @cached_property
    def _get_n_nodal_dofs(self):
        return self.parent_fets.n_nodal_dofs * 6
    
    n_e_dofs = Property(Int, depends_on = 'parent_fets.n_e_dofs')
    @cached_property
    def _get_n_e_dofs(self):
        return self.parent_fets.n_e_dofs * 6
 
               
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
        p_nodal_dofs =  self.parent_fets.n_e_dofs/self.parent_fets.n_nodal_dofs
        p_value = sign(r_ls_value)
        N_e_list = [p_N_mtx[:,i*4:i*4+2]* \
                    (p_value-sign(node_ls_values[i]))\
                 for i in range(0,n_nodes)]
        N_e_mtx = hstack(N_e_list)
        N_enr_mtx = zeros((2,self.parent_fets.n_e_dofs*6), dtype = 'float_')
        for i in range(0,p_nodal_dofs): 
            N_enr_mtx[0,i] = p_N_mtx[0,i]#s
            N_enr_mtx[1,i] = p_N_mtx[0,i]
            
            N_enr_mtx[0,i+2] = p_N_mtx[0,i]#m
            
            N_enr_mtx[1,i+4] = p_N_mtx[0,i]#f
            
            N_enr_mtx[0,i+6] = N_e_mtx[0,i]#sx
            N_enr_mtx[1,i+6] = N_e_mtx[0,i]#sx
            
            N_enr_mtx[0,i+8] = N_e_mtx[0,i]#mx
            
            N_enr_mtx[1,i+10] = N_e_mtx[0,i]#fx
        return N_enr_mtx

    def get_dNr_psi_mtx(self,r_pnt, node_ls_values, r_ls_value):
        '''
        Return the derivatives of the shape functions
        '''
        p_dNr_mtx = self.get_parent_dNr_mtx(r_pnt)
        p_value = sign(r_ls_value)
        n_nodes = self.parent_fets.n_e_dofs/self.parent_fets.n_nodal_dofs
        dNr_mtx = array( [ p_dNr_mtx[:,i]* \
                           (p_value\
                           - sign(node_ls_values[i])) for i in range(0,n_nodes) ])
        return dNr_mtx.T
    
    def get_B_mtx( self, r_pnt, X_mtx, node_ls_values, r_ls_value ):
        J_mtx = self.get_J_mtx(r_pnt,X_mtx)
        dNr_mtx = self.get_parent_dNr_mtx( r_pnt )
        dNr_psi_mtx = self.get_dNr_psi_mtx( r_pnt, node_ls_values, r_ls_value )
        dNx_mtx = dot( inv( J_mtx ), dNr_mtx  )
        dNx_psi_mtx = dot( inv( J_mtx ), dNr_psi_mtx  )
        if dNr_mtx.shape[0] == 1:#1D
            Bx_mtx = zeros( (3,  self.parent_fets.n_e_dofs*6), dtype = 'float_' )
            for i in range(0,(self.parent_fets.n_e_dofs/self.parent_fets.n_nodal_dofs)):
                Bx_mtx[0,i]   = dNx_mtx[0,i]#s  
                Bx_mtx[1,i]   = dNx_mtx[0,i]  
                
                Bx_mtx[0,i+2] = dNx_mtx[0,i]#m
                
                Bx_mtx[1,i+4] = dNx_mtx[0,i]#f
                      
                Bx_mtx[0,i+6] = dNx_psi_mtx[0,i]#sx
                Bx_mtx[1,i+6] = dNx_psi_mtx[0,i]#sx
                
                Bx_mtx[0,i+8] = dNx_psi_mtx[0,i]#mx
                
                Bx_mtx[1,i+10] = dNx_psi_mtx[0,i]#fx      
 
            Bx_mtx[2,:] = [0,0,1,1,-1,-1,0,0,1,1,-1,-1]#TODO:just for 2 nodes
 
        elif dNr_mtx.shape[0] == 2:#2D
            Bx_mtx = zeros( (8,  self.parent_fets.n_e_dofs*5), dtype = 'float_' )
            for i in range(0,(self.parent_fets.n_e_dofs/self.parent_fets.n_nodal_dofs)):
                
                Bx_mtx[0,i*4]   = dNx_mtx[0,i]
                Bx_mtx[1,i*4+1] = dNx_mtx[1,i]
                Bx_mtx[2,i*4]   = dNx_mtx[1,i]
                Bx_mtx[2,i*4+1] = dNx_mtx[0,i]
            
                Bx_mtx[3,i*4+2] = dNx_mtx[0,i]
                Bx_mtx[4,i*4+3] = dNx_mtx[1,i]
                Bx_mtx[5,i*4+2] = dNx_mtx[1,i]
                Bx_mtx[5,i*4+3] = dNx_mtx[0,i]
               
                Bx_mtx[0,i*2+16] = dNx_psi_mtx[0,i]
                Bx_mtx[1,i*2+17] = dNx_psi_mtx[1,i]
                Bx_mtx[2,i*2+16] = dNx_psi_mtx[1,i]
                Bx_mtx[2,i*2+17] = dNx_psi_mtx[0,i]  
            Bx_mtx[6,:] = [1.,0.,-1.,0.,#TODO:just for four nodes
                           1.,0.,-1.,0.,
                           1.,0.,-1.,0., 
                           1.,0.,-1.,0.,
                           1.,0.,-1.,0., 
                           1.,0.,-1.,0.]
            Bx_mtx[7,:] = [0.,1.,0.,-1.,
                           0.,1.,0.,-1.,
                           0.,1.,0.,-1.,
                           0.,1.,0.,-1.,
                           0.,1.,0.,-1.,
                           0.,1.,0.,-1.]    
        return Bx_mtx
