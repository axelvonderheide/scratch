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
from crease_pattern import CreasePattern
from crease_pattern_view import CreasePatternView

def rhombcp():

    cp = CreasePattern()

    cp.nodes = [[ 0, 0, 0 ],
                [ 2, 0, 0 ],
                [ 4, 0, 0 ],
                [ 0, 1, 0 ],
                [ 1, 1, 0 ],
                [ 3, 1, 0 ],
                [ 4, 1, 0]]

    cp.crease_lines = [[ 0, 1 ],
                       [ 1, 2 ],
                       [ 0, 3 ],
                       [ 2, 6 ],
                       [ 0, 4 ],
                       [ 1, 4 ],
                       [ 3, 4 ],
                       [ 5, 4 ],
                       [ 1, 5 ],
                       [ 2, 5 ],
                       [ 6, 5 ],
                        ]

    cp.constraints = [[3, 2],
                      [3, 1],
                      [4, 1],
                      [4, 2],
                      [5, 1],
                      [5, 2],
                      [6, 0],
                      [6, 1],
                      [6, 2],
                      [0, 2]]

    # lift node 0 in z-axes
    cp.constraint_values = np.zeros((10,), dtype = float)
    cp.constraint_values[3] = 0.4
    cp.constraint_values[5] = 0.4
    cp.constraint_values[8] = 0.0
    cp.constraint_values[9] = 0.0

    X = np.zeros((cp.n_dofs,), dtype = float)
    X[3] = 0.1
    X[5] = 0.1

    print 'initial lengths\n', cp.c_lengths
    print 'initial vectors\n', cp.c_vectors

    print 'initial R\n', cp.get_R(X)
    print 'initial dR\n', cp.get_dR(X)

    print 'constrained dofs\n', cp.cnstr_dofs
    print 'free dofs\n', cp.free_dofs

    # Newton-Raphson iteration
    MAX_ITER = 150
    TOLERANCE = 1e-10
    n_steps = 3
    cv = np.copy(cp.constraint_values)

    for k in range(n_steps):
        print 'step', k
        #cp.set_next_node(X)
        cp.constraint_values = (k + 1.) / float(n_steps) * cv
        i = 0
        while i in range(MAX_ITER):
            dR = cp.get_dR(X)[:, cp.free_dofs ]
            print 'dR', dR.shape
            R = cp.get_R(X)
            if np.linalg.norm(R) < TOLERANCE:
                print '==== converged in ', i, 'iterations ===='
                cp.set_next_node(X)
                break
            dX = np.linalg.solve(dR, -R)
            X[ cp.free_dofs ] += dX
            i += 1
        else:
            raise ValueError
            break


    print '========== results =============='
    print 'solution X\n', X
    print 'final positions\n', cp.get_new_nodes(X)
    print 'final vectors\n', cp.get_new_vectors(X)
    print 'final lengths\n', cp.get_new_lengths(X)

    # initialise View
    my_model = CreasePatternView(data = cp)

    my_model.configure_traits()


if __name__ == '__main__':
    rhombcp()
