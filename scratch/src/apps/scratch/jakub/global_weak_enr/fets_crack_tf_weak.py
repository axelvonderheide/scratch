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

class FETSCrackTFWeak( FETSLSEval ):

    x_slice = slice( 0, 4 )

    n_nodal_dofs = Property( Int, depends_on = 'parent_fets.n_nodal_dofs' )
    @cached_property
    def _get_n_nodal_dofs( self ):
        return self.parent_fets.n_nodal_dofs * 2

    n_e_dofs = Property( Int, depends_on = 'parent_fets.n_e_dofs' )
    @cached_property
    def _get_n_e_dofs( self ):
        return self.parent_fets.n_e_dofs * 2


    #---------------------------------------------------------------------
    # Method required to represent the element geometry
    #---------------------------------------------------------------------

    #---------------------------------------------------------------------
    # Method delivering the shape functions for the field variables and their derivatives
    #---------------------------------------------------------------------

    def get_N_mtx( self, r_pnt, node_ls_values, r_ls_value ):
        '''
        Returns the matrix of the shape functions used for the field approximation
        containing zero entries. The number of rows corresponds to the number of nodal
        dofs. The matrix is evaluated for the specified local coordinate r.
        '''
        p_N_mtx = self.parent_fets.get_N_mtx( r_pnt )
        n_nodes = self.parent_fets.n_e_dofs / self.parent_fets.n_nodal_dofs
        p_nodal_dofs = 4
        p_value = sign( r_ls_value )
        psi_weak = self._get_Psi( node_ls_values, p_N_mtx )
        N_e_mtx = hstack( [hstack( ( p_N_mtx[:, i * 4:i * 4 + 2] * \
                          ( p_value - sign( node_ls_values[i] ) ), \
                          p_N_mtx[:, i * 4 + 2:i * 4 + 4] * \
                           psi_weak ) )\
                          for i in range( 0, n_nodes )] )
        N_enr_mtx = hstack( ( p_N_mtx, N_e_mtx ) )
        return N_enr_mtx

    def get_dNr_psi_mtx( self, r_pnt, node_ls_values, r_ls_value ):
        '''
        Return the derivatives of the shape functions
        '''
        p_dNr_mtx = self.get_dNr_mtx( r_pnt )
        p_N_mtx = self.parent_fets.get_N_mtx( r_pnt )
        p_N_red = vstack( ( p_N_mtx[0, ::4], p_N_mtx[0, ::4] ) )#hack for this specific case

        p_value = sign( r_ls_value )
        n_nodes = self.parent_fets.n_e_dofs / self.parent_fets.n_nodal_dofs
        dNr_mtx = hstack( [ p_dNr_mtx[:, i:i + 1] * \
                           ( p_value\
                           - sign( node_ls_values[i] ) ) for i in range( 0, n_nodes ) ] )

        first_mtx = p_N_red * self._get_dPsir( node_ls_values, p_N_mtx, p_dNr_mtx )
        second_mtx = p_dNr_mtx * self._get_Psi( node_ls_values, p_N_mtx )#??
        dNr_weak_mtx = first_mtx + second_mtx
        return hstack( ( dNr_mtx, dNr_weak_mtx ) )

    def _get_Psi( self, node_ls_values, p_N_mtx ):
        p_nodal_dofs = self.parent_fets.n_nodal_dofs
        first = ( hstack( [node_ls_values[i] * \
                   p_N_mtx[:, i * p_nodal_dofs:i * p_nodal_dofs + p_nodal_dofs] \
                   for i in range( 0, self.n_nodes )] ) ).sum()

        third = ( hstack( [p_N_mtx[:, i * p_nodal_dofs:i * p_nodal_dofs + p_nodal_dofs] * \
                     abs( node_ls_values[i] )\
                     for i in range( 0, self.n_nodes )] ) ).sum()
        #print "Psi ", third - abs(first)
        return third - abs( first )

    def _get_dPsir( self, node_ls_values, p_N_mtx, p_dNr_mtx ):
        p_N_red = vstack( ( p_N_mtx[0, ::4], p_N_mtx[0, ::4] ) )#hack for this specific case

        first = sum( hstack( [node_ls_values[i] * \
                   p_N_red[:, i] \
                   for i in range( 0, self.n_nodes ) ] ) )
        second = sum( hstack( [abs( node_ls_values[i] ) * \
                    p_dNr_mtx[:, i] \
                    for i in range( 0, self.n_nodes ) ] ) )
        fourth = sum( hstack( [ p_dNr_mtx[:, i] * \
                    node_ls_values[i]\
                    for i in range( 0, self.n_nodes ) ] ) )
        #print "dPsi ", second - sign(first)*fourth
        return second - sign( first ) * fourth


    def get_B_mtx( self, r_pnt, X_mtx, node_ls_values, r_ls_value ):
        J_mtx = self.get_J_mtx( r_pnt, X_mtx )
        dNr_mtx = self.get_dNr_mtx( r_pnt )
        dNr_psi_mtx = self.get_dNr_psi_mtx( r_pnt, node_ls_values, r_ls_value )
        dNx_mtx = dot( inv( J_mtx ), dNr_mtx )
        dNx_psi_mtx = dot( inv( J_mtx ), dNr_psi_mtx )
        n_nodes = self.parent_fets.n_e_dofs / self.parent_fets.n_nodal_dofs
        if dNr_mtx.shape[0] == 1:#1D#TODO:
            Bx_mtx = zeros( ( 1, self.n_e_dofs ), dtype = 'float_' )
            for i in range( 0, n_nodes ):
                Bx_mtx[0, i] = dNx_mtx[0, i]
                Bx_mtx[0, i + 2] = dNx_psi_mtx[0, i]

        elif dNr_mtx.shape[0] == 2:#2D
            Bx_mtx = zeros( ( 8, self.n_e_dofs ), dtype = 'float_' )
            for i in range( 0, n_nodes ):

                Bx_mtx[0, i * 4] = dNx_mtx[0, i]
                Bx_mtx[1, i * 4 + 1] = dNx_mtx[1, i]
                Bx_mtx[2, i * 4] = dNx_mtx[1, i]
                Bx_mtx[2, i * 4 + 1] = dNx_mtx[0, i]

                Bx_mtx[3, i * 4 + 2] = dNx_mtx[0, i]
                Bx_mtx[4, i * 4 + 3] = dNx_mtx[1, i]
                Bx_mtx[5, i * 4 + 2] = dNx_mtx[1, i]
                Bx_mtx[5, i * 4 + 3] = dNx_mtx[0, i]

                Bx_mtx[0, i * 2 + self.parent_fets.n_e_dofs] = dNx_psi_mtx[0, i]
                Bx_mtx[1, i * 2 + self.parent_fets.n_e_dofs + 1] = dNx_psi_mtx[1, i]
                Bx_mtx[2, i * 2 + self.parent_fets.n_e_dofs] = dNx_psi_mtx[1, i]
                Bx_mtx[2, i * 2 + self.parent_fets.n_e_dofs + 1] = dNx_psi_mtx[0, i]

                Bx_mtx[6, i * 4:i * 4 + 4] = [1., 0., -1., 0.]#TODO:this is just 2d
                Bx_mtx[7, i * 4:i * 4 + 4] = [0., 1., 0., -1.]
        return Bx_mtx
