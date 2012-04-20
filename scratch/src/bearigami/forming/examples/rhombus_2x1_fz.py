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

from enthought.mayavi.core.api import PipelineBase
from enthought.mayavi.core.ui.api import MayaviScene, SceneEditor, \
    MlabSceneModel
from enthought.mayavi.modules.axes import Axes

from enthought.traits.api import HasTraits, Range, Instance, on_trait_change, \
    Trait, Property, Constant, DelegatesTo, cached_property, Str, Delegate, \
    Button, Int
from enthought.traits.ui.api import View, Item, Group, ButtonEditor
from enthought.mayavi import mlab
from enthought.mayavi.core.api import Engine
import numpy as np

# own Modules
from bearigami.forming.crease_pattern import CreasePattern
from bearigami.forming.crease_pattern_view import CreasePatternView

def run():

    cp = CreasePattern(n_steps = 30, MAX_ITER = 1000)

    cp.nodes = [[ 0, 0, 0 ],
                [ 2, 0, 0 ],
                [ 4, 0, 0 ],
                [ 6, 0, 0 ],
                [ 8, 0, 0 ],
                [ 0, 1, 0 ],
                [ 1, 1, 0 ],
                [ 3, 1, 0 ],
                [ 5, 1, 0 ],
                [ 7, 1, 0 ],
                [ 8, 1, 0 ],
                ]

    cp.crease_lines = [[ 0, 1 ], #1
                       [ 1, 2 ], #2
                       [ 2, 3 ], #3
                       [ 3, 4 ], #4
                       [ 0, 5 ],
                       [ 0, 6 ],
                       [ 1, 6 ],
                       [ 1, 7 ],
                       [ 2, 7 ],
                       [ 2, 8 ],
                       [ 3, 8 ],
                       [ 3, 9 ],
                       [ 4, 9 ],
                       [ 4, 10],
                       [ 5, 6 ],
                       [ 6, 7 ],
                       [ 7, 8 ],
                       [ 8, 9 ],
                       [ 9, 10 ],
                        ]

    cp.cnstr_lhs = [[(0, 2, 1.0)], # n0.z = 0.0
                    [(4, 2, 1.0)], # n4.z = 0.0
                    [(2, 2, 1.0)], # n2.z = control node
                    [(2, 0, 1.0)], # n2.x = 0
                    [(5, 1, 1.0)], # n5.y = 0
                    [(5, 1, 1.0), (6, 1, -1.0)], # n6.y = n3.y
                    [(5, 1, 1.0), (7, 1, -1.0)], # n7.y = n3.y
                    [(5, 1, 1.0), (8, 1, -1.0)], # n8.y = n3.y
                    [(5, 1, 1.0), (9, 1, -1.0)], # n9.y = n3.y
                    [(5, 1, 1.0), (10, 1, -1.0)], # n10.y = n3.y
                    [(0, 1, 1.0), (1, 1, -1.0)], # n1.y = n0.y
                    [(0, 1, 1.0), (2, 1, -1.0)], # n2.y = n0.y
                    [(0, 1, 1.0), (3, 1, -1.0)], # n3.y = n0.y
                    [(0, 1, 1.0), (4, 1, -1.0)], # n4.y = n0.y
                    ]

    # lift node 0 in z-axes
    cp.cnstr_rhs = np.zeros((14,), dtype = float)
    cp.cnstr_rhs[2] = 3.0

    X = np.zeros((cp.n_dofs,), dtype = float)
    X[2] = 0.0
    X[5] = 0.05
    X[8] = 0.0
    X[11] = 0.05
    X[14] = 0.0
    X[17] = -0.05
    X[20] = 0.05
    X[23] = 0.01
    X[26] = 0.01
    X[29] = 0.05
    X[32] = -0.05

    print 'initial lengths\n', cp.c_lengths
    print 'initial vectors\n', cp.c_vectors

    print 'initial R\n', cp.get_R(X)
    dR = cp.get_dR(X)
    np.savetxt('dR.arr', dR, '%g')
    np.savetxt
    print 'initial dR\n', cp.get_dR(X)

    X = cp.solve(X)

    print '========== results =============='
    print 'solution X\n', X
    print 'final positions\n', cp.get_new_nodes(X)
    print 'final vectors\n', cp.get_new_vectors(X)
    print 'final lengths\n', cp.get_new_lengths(X)

    # initialise View
    my_model = CreasePatternView(data = cp)

    my_model.configure_traits()


if __name__ == '__main__':
    run()
