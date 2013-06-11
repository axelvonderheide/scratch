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
from bearigami.forming.rhombus_crease_pattern import RhombusCreasePattern
from bearigami.forming.crease_pattern_view import CreasePatternView


def cp01(L_x = 4, L_y = 4, n_x = 1, n_y = 2,
         n_steps = 100, skew_coeff = 0.0):
    '''Exploit symmetric constraints
    '''
    cp = RhombusCreasePattern(n_steps = n_steps,
                              L_x = L_x,
                              L_y = L_y,
                              n_x = n_x,
                              n_y = n_y,
                              show_iter = False,
                              MAX_ITER = 500)

    n_h = cp.n_h
    n_i = cp.n_i
    n_v = cp.n_v
    n_h_idx = n_y / 4

    x_links = []
    y_links = []
    z_links = []

    z_nodes = n_h[(0, 0, -1, -1), (0, -1, -1, 0)].flatten()
    print 'z_nodes', z_nodes

    #z_cnstr = [[(n, 2, 1.0)] for n in z_nodes]
    x_cnstr = [[(n_h[0, 0], 0, 1.0)]]
    y_cnstr = [[(n_h[0, 0], 1, 1.0)]]
    z_cnstr = [[(n_h[0, 0], 2, 1.0)]]

    for n_arr in n_h[:, (0, -1)].T:
        for idx, n in enumerate(n_arr[1:]):
            n_x = len(n_arr)
            coeff = skew_coeff * float(idx + 1) / float(n_x)
            y_links.append([(n_arr[0], 1, 1.0 - coeff), (n, 1, -1.0)])

    for n in n_h[0, 1:]:
        z_links.append([(n_h[0, 0], 2, 1.0), (n, 2, -1.0)])
        x_links.append([(n_h[0, 0], 0, 1.0), (n, 0, -1.0)])
    #x_links.append([(n_h[0, -1], 0, 1.0), (n_h[0, -1], 1, -0.5)])

    for n in n_v[-1, 1:]:
        x_links.append([(n_v[-1, 0], 0, 1.0), (n, 0, -1.0)])

    for n0, n1 in zip(n_v[0, :], n_v[-1, :]):
        z_links.append([(n0, 2, 1.0), (n1, 2, -1.0)])

    #cntrl = [[(n_h[-1, -1], 1, 1.0)]]
    cntrl = [[(n_h[-1, 1], 0, 1.0)]]

    print 'x_cnstr', len(x_cnstr)
    print 'y_cnstr', len(y_cnstr)
    print 'z_cnstr', len(z_cnstr)
    print 'x_links', len(x_links)
    print 'y_links', len(y_links)
    print 'z_links', len(z_links)

    cp.cnstr_lhs = z_cnstr + x_links + y_links + z_links + x_cnstr + y_cnstr + cntrl
    #cp.cnstr_lhs = z_cnstr

    # lift node 0 in z-axes
    cp.cnstr_rhs = np.zeros((len(cp.cnstr_lhs),), dtype = float)
    cp.cnstr_rhs[-1] = -L_x * 0.34

#    cp.cnstr_rhs[-1] = -L_y * 0.9999

    return cp

if __name__ == '__main__':

    cp = cp01(L_x = 14, L_y = 8, n_x = 4, n_y = 4,
              n_steps = 40, skew_coeff = 0.0)
    X0 = cp.generate_X0()
    #cp.set_next_node(X0)

    print 'n_dofs', cp.n_dofs
    print 'n_c', cp.n_c
    print 'necessary constraints', cp.n_dofs - cp.n_c
    print 'cnstr', len(cp.cnstr_lhs)

    X_vct = cp.solve(X0)

#    print 'new nodes'
#    print cp.get_new_nodes(X_vct)
#    print 'new lengths'
#    print cp.get_new_lengths(X_vct)

    # initialise View
    my_model = CreasePatternView(data = cp, show_cnstr = True)
    my_model.configure_traits()

