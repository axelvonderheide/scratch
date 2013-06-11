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
# Created on Sep 8, 2011 by: matthias


from enthought.traits.api import HasTraits, Property, cached_property, Event, \
    Array, Instance, Int, Directory, Range, on_trait_change, Bool, Trait, Constant, \
    Tuple, Interface, implements, Enum

from enthought.traits.trait_types import DelegatesTo
from enthought.traits.ui.api import Item, View, HGroup, RangeEditor

import numpy as np
from crease_pattern import CreasePattern

class CreasePatternTest(HasTraits):



    def cp_t1(self):

        #------------------------------------------------------------------------------ 
        # test with 2 creaselines
        #------------------------------------------------------------------------------ 

        cp = CreasePattern()

        cp.nodes = [[ 0, 0, 0 ],
                    [ 1, 0, 0 ],
                    [ 1, 1, 0 ]]

        cp.crease_lines = [[ 0, 1 ], [1, 2]] # , [2, 0]]

        cp.constraints = [[0, 0], [0, 1], [0, 2],
                          [2, 0], [2, 1], [2, 2],
                          [1, 2]
                          ]

        cp.constraint_values = [0, 0, 0, 0, 0, 0, 0.5]

        X = np.zeros((cp.n_dofs,), dtype = float)

        print 'initial lengths\n', cp.c_lengths
        print 'initial vectors\n', cp.c_vectors

        print 'initial R\n', cp.get_R(X)
        print 'initial dR\n', cp.get_dR(X)

        print 'constrained dofs\n', cp.cnstr_dofs
        print 'free dofs\n', cp.free_dofs

        # Newton-Raphson iteration
        MAX_ITER = 150
        TOLERANCE = 1e-10

        for i in range(MAX_ITER):
            dR = cp.get_dR(X)[:, cp.free_dofs ]
            R = cp.get_R(X)
            if np.linalg.norm(R) < TOLERANCE:
                print '==== converged in ', i, 'iterations ===='
                break
            dX = np.linalg.solve(dR, -R)
            X[ cp.free_dofs ] += dX

        print '========== results =============='
        print 'solution X\n', X
        print 'final positions\n', cp.get_new_nodes(X)
        print 'final vectors\n', cp.get_new_vectors(X)
        print 'final lengths\n', cp.get_new_lengths(X)

    def cp_t2(self):

        #------------------------------------------------------------------------------ 
        # test with 3 creaselines
        #------------------------------------------------------------------------------ 

        cp = CreasePattern()

        cp.nodes = [[ 0, 0, 0 ],
                    [ 1, 0, 0 ],
                    [ 1, 1, 0 ]]

        cp.crease_lines = [[ 0, 1 ], [1, 2], [2, 0]]

        cp.constraints = [[0, 0], [0, 1], [0, 2],
                          [1, 2], [2, 2], [1, 1],
                          ]

        cp.constraint_values = [0, 0, 0, 0, 0, 0.1]

        X = np.zeros((cp.n_dofs,), dtype = float)

        print 'initial lengths\n', cp.c_lengths
        print 'initial vectors\n', cp.c_vectors

        print 'initial R\n', cp.get_R(X)
        print 'initial dR\n', cp.get_dR(X)

        print 'constrained dofs\n', cp.cnstr_dofs
        print 'free dofs\n', cp.free_dofs

        # Newton-Raphson iteration
        MAX_ITER = 150
        TOLERANCE = 1e-10

        for i in range(MAX_ITER):
            dR = cp.get_dR(X)[:, cp.free_dofs ]
            R = cp.get_R(X)
            if np.linalg.norm(R) < TOLERANCE:
                print '==== converged in ', i, 'iterations ===='
                break
            dX = np.linalg.solve(dR, -R)
            X[ cp.free_dofs ] += dX

        print '========== results =============='
        print 'solution X\n', X
        print 'final positions\n', cp.get_new_nodes(X)
        print 'final vectors\n', cp.get_new_vectors(X)
        print 'final lengths\n', cp.get_new_lengths(X)


    def cp_t3(self):

        #------------------------------------------------------------------------------ 
        # test 6 triangles (standard pattern)
        # 7 nodes                              
        # 12 creaselines
        # 3 * 7 - 12 constrains needed
        #------------------------------------------------------------------------------ 

        cp = CreasePattern()

        cp.nodes = [[ 0, 0, 0 ],
                    [ 1, 0, 0 ],
                    [ 1, 1, 0 ],
                    [ -1, 1, 0 ],
                    [ -1, 0, 0 ],
                    [ -1, -1, 0 ],
                    [ 1, -1, 0 ]]

        cp.crease_lines = [[ 0, 1 ],
                           [ 1, 2 ],
                           [ 2, 0 ],
                           [ 2, 3 ],
                           [ 3, 0 ],
                           [ 3, 4 ],
                           [ 4, 0 ],
                           [ 4, 5 ],
                           [ 5, 0 ],
                           [ 5, 6 ],
                           [ 6, 0 ],
                           [ 6, 1 ]]


        cp.constraints = [[0, 1],
                          [0, 2],
                          [1, 0],
                          [1, 1],
                          [1, 2],
                          [4, 1],
                          [2, 0],
                          [6, 0],
                          [4, 2],
                          ]

        # lift node 0 in z-axes
        cp.constraint_values = [0, 0.5, 0, 0, 0, 0, 0, 0, 0]

        X = np.zeros((cp.n_dofs,), dtype = float)

        print 'initial lengths\n', cp.c_lengths
        print 'initial vectors\n', cp.c_vectors

        print 'initial R\n', cp.get_R(X)
        print 'initial dR\n', cp.get_dR(X)

        print 'constrained dofs\n', cp.cnstr_dofs
        print 'free dofs\n', cp.free_dofs

        # Newton-Raphson iteration
        MAX_ITER = 150
        TOLERANCE = 1e-10
        n_steps = 8
        cv = np.copy(cp.constraint_values)

        for k in range(n_steps):
            print 'step', k
            cp.constraint_values = (k + 1.) / float(n_steps) * cv
            for i in range(MAX_ITER):
                dR = cp.get_dR(X)[:, cp.free_dofs ]
                print 'dR', dR.shape
                R = cp.get_R(X)
                if np.linalg.norm(R) < TOLERANCE:
                    print '==== converged in ', i, 'iterations ===='
                    break
                dX = np.linalg.solve(dR, -R)
                X[ cp.free_dofs ] += dX
            cp.set_next_node(X)

        print '========== results =============='
        print 'solution X\n', X
        print 'final positions\n', cp.get_new_nodes(X)
        print 'final vectors\n', cp.get_new_vectors(X)
        print 'final lengths\n', cp.get_new_lengths(X)

        # initialise View
        from crease_pattern_view import CreasePatternView
        my_model = CreasePatternView(data = cp)
        print my_model.data.iteration_nodes.shape

        my_model.configure_traits()


if __name__ == '__main__':
    cpt = CreasePatternTest()

    cpt.cp_t3()
