#-------------------------------------------------------------------------------
#
# Copyright (c) 2009, IMB, RWTH Aachen.
# All rights reserved.
#
# This software is provided without warranty under the terms of the BSD
# license included in simvisage/LICENSE.txt and may be redistributed only
# under the conditions described in the aforementioned license.  The license
# is also available online at http://www.simvisage.com/licenses/BSD.txt
#
# Thanks for using Simvisage open source!
#
# Created on Sep 7, 2011 by: rch

'''
Created on 25.08.2011

@author: schmerl
'''
from enthought.traits.api import HasTraits, Property, cached_property, Event, \
    Array, Instance, Int, Directory, Range, on_trait_change, Bool, Trait, Constant, \
    Tuple, Interface, implements, Enum, List, Float

from enthought.traits.trait_types import DelegatesTo
from enthought.traits.ui.api import Item, View, HGroup, RangeEditor

import numpy as np

from bearigami.forming.cnstr_entity import Cnstr_Entity, Cnstr_Line, Cnstr_Plane

class CreasePattern(HasTraits):
    '''
        Structure of triangulated Crease-Patterns
    '''
    #===============================================================================
    # Input data structure 
    #===============================================================================

    # all nodes in X,Y,Z array    
    nodes = Array(value = [], dtype = float)

    # all crease lines as index-table
    crease_lines = Array(value = [], dtype = int)

    # constrained node indices
    # define the pairs (node, dimension) affected by the constraint
    # stored in the constrained_x array
    #
    # define the constraint in the form
    # crrs_lhs = [ [(node_1, dir_1, coeff_1),(node_2, dir_2, coeff_2)], # first constraint
    #              [(node_i, dir_i, coeff_i)], # second constraint
    # ctrs_rhs = [ value_first, velue_second ]
    # 
    # left-hand side coefficients of the constraint equations 
    cnstr_lhs = List()
    # right-hand side values of the constraint equations
    cnstr_rhs = Array(value = [], dtype = float)
    # list of Constrain-Objects
    cnstr = Array(value = [])

    #===============================================================================
    # Enumeration of dofs 
    #===============================================================================

    all_dofs = Property(Array, depends_on = 'constraints')
    @cached_property
    def _get_all_dofs(self):
        return np.arange(self.n_dofs).reshape(self.n_n, self.n_d)

    #===============================================================================
    # Convenience properties providing information about the input 
    #===============================================================================
    n_n = Property
    def _get_n_n(self):
        '''Number of crease nodes'''
        return self.nodes.shape[0]

    n_c = Property
    def _get_n_c(self):
        '''Number of crease lines'''
        return self.crease_lines.shape[0]

    n_d = Constant(3)

    # total number of dofs
    n_dofs = Property(depends_on = 'n_d,n_c,n_d')
    @cached_property
    def _get_n_dofs(self):
        return self.n_n * self.n_d

    #===========================================================================
    # Dependent interim results
    #===========================================================================
    c_vectors = Property(Array, depends_on = 'nodes, crease_lines')
    @cached_property
    def _get_c_vectors(self):
        '''
            Calculates the c of the crease lines.
        '''
        n = self.nodes[...]

        cl = self.crease_lines
        return n[ cl[:, 1] ] - n[ cl[:, 0] ]

    c_lengths = Property(Array, depends_on = 'nodes, crease_lines')
    @cached_property
    def _get_c_lengths(self):
        '''
            Calculates the lengths of the crease lines.
        '''
        c = self.c_vectors
        return np.sqrt(np.sum(c ** 2, axis = 1))

    #===============================================================================
    # Verification procedures to check the compliance with the constant length criteria. 
    #===============================================================================
    def get_new_nodes(self, X_vct):
        '''
            Calculates the lengths of the crease lines.
        '''
        X = X_vct.reshape(self.n_n, self.n_d)
        return self.nodes + X

    def get_new_vectors(self, X_vct):
        '''
            Calculates the lengths of the crease lines.
        '''
        cX = self.get_new_nodes(X_vct)
        cl = self.crease_lines
        return cX[ cl[:, 1] ] - cX[ cl[:, 0] ]

    def get_new_lengths(self, X_vct):
        '''
            Calculates the lengths of the crease lines.
        '''
        cV = self.get_new_vectors(X_vct)
        return np.sqrt(np.sum(cV ** 2, axis = 1))

    #===============================================================================
    # Get the predictor and corrector.
    #===============================================================================

    def get_length_R(self, X_vct):
        ''' Calculate the residuum for constant crease length
        given the fold vector dX.
        '''
        j = self.crease_lines[:, 1]
        i = self.crease_lines[:, 0]

#        X_vct[ self.cnstr_dofs ] = self.constraint_values
        X = X_vct.reshape(self.n_n, self.n_d)
        Xj = X[j]
        Xi = X[i]

        CXi = np.sum(self.c_vectors * Xi, axis = 1)
        CXj = np.sum(self.c_vectors * Xj, axis = 1)

        Xij = np.sum(Xi * Xj, axis = 1)
        Xii = np.sum(Xi ** 2, axis = 1)
        Xjj = np.sum(Xj ** 2, axis = 1)
        R = 2 * CXj - 2 * CXi - 2 * Xij + Xii + Xjj

        return R

    def get_length_dR(self, X_vct):
        ''' Calculate the jacobian of the residuum at the instantaneous
        configuration dR
        '''
        i = self.crease_lines[:, 0]
        j = self.crease_lines[:, 1]

        X = X_vct.reshape(self.n_n, self.n_d)
        Xj = X[j]
        Xi = X[i]

        dR_i = -2 * self.c_vectors + 2 * Xi - 2 * Xj
        dR_j = 2 * self.c_vectors + 2 * Xj - 2 * Xi

        dR = np.zeros((self.n_c, self.n_n, self.n_d), dtype = float)

        # running crease line index
        cidx = np.arange(self.n_c)
        dR[ cidx, i, : ] += dR_i
        dR[ cidx, j, : ] += dR_j

        # reshape the 3D matrix to a 2D matrix 
        # with rows for crease lines and columns representing 
        # the derivatives with respect to the node displacements
        # in 3d.
        # 
        return dR.reshape(self.n_c, self.n_n * self.n_d)

    def get_cnstr_R(self, X_vct):
        ''' Calculate the residuum for given constraint equations
        '''
        X = X_vct.reshape(self.n_n, self.n_d)
        Rc = np.zeros((len(self.cnstr_lhs),))
        for i, cnstr in enumerate(self.cnstr_lhs):
            for n, d, c in cnstr:
                Rc[i] += c * X[n, d] - self.cnstr_rhs[i]
        return Rc

    def get_cnstr_dR(self, X_vct):
        ''' Calculate the residuum for given constraint equations
        '''
        dRc = np.zeros((len(self.cnstr_lhs), self.n_dofs))
        for i, cnstr in enumerate(self.cnstr_lhs):
            for n, d, c in cnstr:
                dof = 3 * n + d
                dRc[i, dof] += c
        return dRc

    def get_R(self, X_vct):
        return np.hstack([self.get_length_R(X_vct),
                          self.get_cnstr_R(X_vct),
                          ])

    def get_dR(self, X_vct):
        return np.vstack([self.get_length_dR(X_vct),
                          self.get_cnstr_dR(X_vct),
                          ])

    #===========================================================================
    # Folding algorithm - Newton-Raphson
    #===========================================================================

    MAX_ITER = Int(50)
    TOLERANCE = Float(1e-10)
    n_steps = Int(20)
    show_iter = Bool(False)

    def solve(self, X):

        # Newton-Raphson iteration
        MAX_ITER = self.MAX_ITER
        TOLERANCE = self.TOLERANCE
        n_steps = self.n_steps

        cnstr_rhs = np.copy(self.cnstr_rhs)

        for k in range(n_steps):
            print 'step', k,
            #self.set_next_node(X)
            i = 0
            self.cnstr_rhs = (k + 1.) / float(n_steps) * cnstr_rhs

            while i <= MAX_ITER:
                dR = self.get_dR(X)
                R = self.get_R(X)
                nR = np.linalg.norm(R)
                if nR < TOLERANCE:
                    print '==== converged in ', i, 'iterations ===='
                    self.set_next_node(X)
                    break
                dX = np.linalg.solve(dR, -R)
                X += dX
                if self.show_iter:
                    self.set_next_node(X)
                i += 1
            else:
                print '==== did not converge in %d interations ====' % i
                return X

        return X

    #===============================================================================
    # calculate Nodeposition for each iterationstep for visualation
    #===============================================================================

    iteration_nodes = Array(value = [], dtype = float)

    def set_next_node(self, X_vct):
        '''
           Calculates the position of nodes for this iteration.
        '''
        if(self.iteration_nodes.shape == (0,)):
            self.iteration_nodes = [self.nodes]

        X = X_vct.reshape(self.n_n, self.n_d)
        nextnode = self.nodes + X
        self.iteration_nodes = np.vstack((self.iteration_nodes, [nextnode]))


    def get_cnstr_pos(self, nodes = None):
        '''Get the coordinates of the constraints.
        '''
        if nodes == None:
            nodes = self.nodes

        cnstr_pos_list = []
        for lhs, rhs in zip(self.cnstr_lhs, self.cnstr_rhs):
            # get the node index
            for n_idx, dir, coeff in lhs:
                X_n = nodes[n_idx, :]
                U_n = np.zeros((self.n_d,), dtype = 'f')
                U_n[ dir ] = 1.0
                cnstr_pos = np.hstack([X_n, U_n])
                cnstr_pos_list.append(cnstr_pos)
        return np.vstack(cnstr_pos_list)
    
    

if __name__ == '__main__':

    # trivial example with a single truss positioned 
    # alone x-axis in the range [0,1] 
    # the first node can moves along the the y-axis
    # the second node can move along the line y = 1 - x
    # 
    cp = CreasePattern()

    cp.nodes = [[ 0, 0, 0 ],
                [ 1, 0, 0 ]]

    cp.crease_lines = [[ 0, 1 ]]

    cp.cnstr_lhs = [
                    [(0, 0, 1.0)],
                    [(0, 2, 1.0)],
                    [(1, 0, 1.0)],
                    [(1, 0, 1.0), (1, 1, 1.0)],
                    [(1, 2, 1.0)]
                    ]

    cp.cnstr_rhs = [0.0, 0.0, -0.1, 0.0, 0.0]

    X = np.zeros((cp.n_dofs,), dtype = float)

    # NOTE: there must be a nonzero initial value
    # of the y-coordinate of one of the nodes
    # otherwise the tangential matrix is singular.
    #  
    X[4] += 0.001

    print 'initial lengths\n', cp.c_lengths
    print 'initial vectors\n', cp.c_vectors

    print 'initial R\n', cp.get_length_R(X)
    print 'initial dR\n', cp.get_length_dR(X)

    print 'cnstr Rc\n', cp.get_cnstr_R(X)
    print 'cnstr dRc\n', cp.get_cnstr_dR(X)

    print 'R\n', cp.get_R(X)
    print 'dR\n', cp.get_dR(X)

    print 'constraint positions', cp.get_cnstr_pos()

    X = cp.solve(X)

    print '========== results =============='
    print 'solution X\n', X
    print 'final positions\n', cp.get_new_nodes(X)
#    print 'final vectors\n', cp.get_new_vectors(X)
    print 'final lengths\n', cp.get_new_lengths(X)
