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
# Created on May 11, 2011 by: rch

from numpy import array, vstack, zeros
from enthought.tvtk.api import tvtk
from enthought.mayavi import mlab

if __name__ == '__main__':

    points = array( [[ 1., 1.   , 0.   ],
                     [-1., 1.   , 0.   ],
                     [-1., -0.31 , 0.   ],
                     [ 1., -0.31 , 0.
                     ]], dtype = 'float_' )

    print points[ -1:, :]
    p_a = vstack( [points[-1:, :], points[:-1, :]] )
    p_b = vstack( [points[:1, :], points[1:, :] ] )

    points_m = ( p_a + p_b ) / 2
    print points_m

    points_x = zeros( ( points.shape[0] * 2, points.shape[1] ), dtype = 'float' )
    points_x[::2, :] = points[:, :]
    points_x[1::2, :] = points_m[:, :]

    print points_x

    # Create a polydata with the points we just created.
    profile = tvtk.PolyData( points = points_x )

    # Perform a 2D Delaunay triangulation on them.
    delny = tvtk.Delaunay2D( input = profile, offset = 1.e1 )

    tri = delny.output
    tri.update()#initiate triangulation
    triangles = array( tri.polys.data, dtype = 'int_' )

    print triangles

    # Author: Gael Varoquaux <gael.varoquaux@normalesup.org>
    # Copyright (c) 2009, Enthought, Inc.
    # License: BSD Style.
    import numpy as np
    # Create data with x and y random in the [-2, 2] segment, and z a
    # Gaussian function of x and y.
    np.random.seed( 12345 )
    x = 4 * ( np.random.random( 500 ) - 0.5 )
    y = 4 * ( np.random.random( 500 ) - 0.5 )
    def f( x, y ):
        return np.exp( -( x ** 2 + y ** 2 ) )
    z = f( x, y )

    mlab.figure( 1, fgcolor = ( 0, 0, 0 ), bgcolor = ( 1, 1, 1 ) )
    # Visualize the points
    pts = mlab.points3d( x, y, z ) # , z, scale_mode = 'none', scale_factor = 0.2 )
    # Create and visualize the mesh
    mesh = mlab.pipeline.delaunay2d( pts )
    surf = mlab.pipeline.surface( mesh )
    mlab.view( 47, 57, 8.2, ( 0.1, 0.15, 0.14 ) )
    mlab.show()
