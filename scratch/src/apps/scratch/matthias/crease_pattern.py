'''
Created on 25.08.2011

@author: schmerl
'''
from enthought.traits.api import HasTraits, Property, cached_property, Event, \
    Array, Instance, Int, Directory, Range, on_trait_change, Bool, Trait, Constant, \
    Tuple, Interface, implements, Enum

from enthought.traits.trait_types import DelegatesTo
from enthought.traits.ui.api import Item, View, HGroup, RangeEditor

import numpy as np

class CreasePattern( HasTraits ):
    '''
        Structure of triangulated Crease-Patterns
    '''
    #===============================================================================
    # Input data structure 
    #===============================================================================

    # all nodes in X,Y,Z array    
    nodes = Array( value = [], dtype = float )

    # all crease lines as index-table
    crease_lines = Array( value = [], dtype = int )

    # constrained node indices
    constraints = Array( value = [], dtype = int )
    # positions of constrained nodes
    constrained_x = Array( value = [], dtype = float )

    #===============================================================================
    # Enumeration of dofs 
    #===============================================================================

    all_dofs = Property( Array, depends_on = 'constraints' )
    @cached_property
    def _get_all_dofs( self ):
        return np.arange( self.n_dofs ).reshape( self.n_n, self.n_d )

    cnstr_dofs = Property( Array, depends_on = 'constraints' )
    @cached_property
    def _get_cnstr_dofs( self ):
        return self.all_dofs[ tuple( self.constraints.T ) ]

    # get the free dofs
    free_dofs = Property( Array, depends_on = 'constraints' )
    @cached_property
    def _get_free_dofs( self ):
        mask = np.repeat( True, self.n_dofs ).reshape( self.n_n, self.n_d )
        mask[ tuple( self.constraints.T ) ] = False
        return self.all_dofs[ mask ]

    #===============================================================================
    # Convenience properties providing information about the input 
    #===============================================================================
    n_n = Property
    def _get_n_n( self ):
        '''Number of crease nodes'''
        return self.nodes.shape[0]

    n_c = Property
    def _get_n_c( self ):
        '''Number of crease lines'''
        return self.crease_lines.shape[0]

    n_d = Constant( 3 )

    # total number of dofs
    n_dofs = Property( depends_on = 'n_d,n_c,n_d' )
    @cached_property
    def _get_n_dofs( self ):
        return self.n_n * self.n_d


    #===========================================================================
    # Dependent interim results
    #===========================================================================
    c_vectors = Property( Array, depends_on = 'nodes, crease_lines' )
    @cached_property
    def _get_c_vectors( self ):
        '''
            Calculates the c of the crease lines.
        '''
        n = self.nodes[...]

        cl = self.crease_lines
        return n[ cl[:, 1] ] - n[ cl[:, 0] ]

    c_lengths = Property( Array, depends_on = 'nodes, crease_lines' )
    @cached_property
    def _get_c_lengths( self ):
        '''
            Calculates the lengths of the crease lines.
        '''
        c = self.c_vectors
        return np.sqrt( np.sum( c ** 2, axis = 1 ) )

    def get_new_nodes( self, X_vct ):
        '''
            Calculates the lengths of the crease lines.
        '''
        X_vct[ self.cnstr_dofs ] = self.constrained_x
        X = X_vct.reshape( self.n_n, self.n_d )
        return self.nodes + X

    def get_new_vectors( self, X_vct ):
        '''
            Calculates the lengths of the crease lines.
        '''
        cX = self.get_new_nodes( X_vct )
        cl = self.crease_lines
        return cX[ cl[:, 1] ] - cX[ cl[:, 0] ]

    def get_new_lengths( self, X_vct ):
        '''
            Calculates the lengths of the crease lines.
        '''
        cV = self.get_new_vectors( X_vct )
        return np.sqrt( np.sum( cV ** 2, axis = 1 ) )

    #===============================================================================
    # Get the predictor and corrector.
    #===============================================================================

    def get_R( self, X_vct ):
        ''' Calculate the residuum for constant crease length
        given the fold vector dX.
        '''
        j = cp.crease_lines[:, 1]
        i = cp.crease_lines[:, 0]

        X_vct[ self.cnstr_dofs ] = self.constrained_x
        X = X_vct.reshape( self.n_n, self.n_d )
        Xj = X[j]
        Xi = X[i]

        CXi = np.sum( cp.c_vectors * Xi, axis = 1 )
        CXj = np.sum( cp.c_vectors * Xj, axis = 1 )

        Xij = np.sum( Xi * Xj, axis = 1 )
        Xii = np.sum( Xi ** 2, axis = 1 )
        Xjj = np.sum( Xj ** 2, axis = 1 )
        R = 2 * CXj - 2 * CXi - 2 * Xij + Xii + Xjj

        return R

    def get_dR( self, X_vct ):
        ''' Calculate the jacobian of the residuum at the instantaneous
        configuration dR
        '''
        i = self.crease_lines[:, 0]
        j = self.crease_lines[:, 1]

        X_vct[ self.cnstr_dofs ] = self.constrained_x
        X = X_vct.reshape( self.n_n, self.n_d )
        Xj = X[j]
        Xi = X[i]

        dR_i = -2 * cp.c_vectors + 2 * Xi - 2 * Xj
        dR_j = 2 * cp.c_vectors + 2 * Xj - 2 * Xi

        dR = np.zeros( ( cp.n_c, cp.n_n, cp.n_d ), dtype = float )

        # running crease line index
        cidx = np.arange( cp.n_c )
        dR[ cidx, i, : ] += dR_i
        dR[ cidx, j, : ] += dR_j

        # reshape the 3D matrix to a 2D matrix 
        # with rows for crease lines and columns representing 
        # the derivatives with respect to the node displacements
        # in 3d.
        # 
        return dR.reshape( cp.n_c, cp.n_n * cp.n_d )

if __name__ == '__main__':

    cp = CreasePattern()

    cp.nodes = [[ 0, 0, 0 ],
                [ 1, 0, 0 ],
                [ 1, 1, 0 ]]

    cp.crease_lines = [[ 0, 1 ], [1, 2]] # , [2, 0]]

    cp.constraints = [[0, 0], [0, 1], [0, 2],
                      [2, 0], [2, 1], [2, 2],
                      [1, 2]
                      ]

    cp.constrained_x = [0, 0, 0, 0, 0, 0, 0.5]

    X = np.zeros( ( cp.n_dofs, ), dtype = float )

    X[4] = -0.3
    X[5] = 0.3

    print 'initial lengths\n', cp.c_lengths
    print 'initial vectors\n', cp.c_vectors

    print 'initial R\n', cp.get_R( X )
    print 'initial dR\n', cp.get_dR( X )

    print 'constrained dofs\n', cp.cnstr_dofs
    print 'free dofs\n', cp.free_dofs

    # Newton-Raphson iteration
    MAX_ITER = 150
    TOLERANCE = 1e-10

    for i in range( MAX_ITER ):
        dR = cp.get_dR( X )[:, cp.free_dofs ]
        R = cp.get_R( X )
        if np.linalg.norm( R ) < TOLERANCE:
            print '==== converged in ', i, 'iterations ===='
            break
        dX = np.linalg.solve( dR, -R )
        X[ cp.free_dofs ] += dX

    print '========== results =============='
    print 'solution X\n', X
    print 'final positions\n', cp.get_new_nodes( X )
    print 'final vectors\n', cp.get_new_vectors( X )
    print 'final lengths\n', cp.get_new_lengths( X )
