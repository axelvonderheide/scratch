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

def rhombcp():

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

    cp.cnstr_rhs = [0.0, 0.0, -1.99, 0.0, 0.0]

    X = np.zeros((cp.n_dofs,), dtype = float)
    X[1] = 0.01

    print 'initial lengths\n', cp.c_lengths
    print 'initial vectors\n', cp.c_vectors

    print 'initial R\n', cp.get_R(X)
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
    rhombcp()
