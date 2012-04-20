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
# Created on Jul 22, 2010 by: rch

from numpy import \
    loadtxt, ones_like, vstack, c_, hstack, array, cumsum, \
    zeros_like, zeros, min as npmin

import wxversion

wxversion.select( '2.8' )

from os.path import join
import os, re

from promod.simdb import SimDB
simdb = SimDB()

data_dir = join( simdb.exdata_dir, 'trc', 'yarn_structure', 'VET', 'raw_data' )

from enthought.tvtk.api import tvtk
from enthought.mayavi.scripts import mayavi2

from enthought.mayavi import mlab


def get_slice_files( datadir ):
    """Given a dirname, returns a list of all slice files."""
    result = []
    group = []
    paths = os.listdir( datadir )  # list of paths in that dir
    for fname in paths:
        match = re.search( r'\w+Schnitt(\d+).txt', fname )
        if match:
            result.append( os.path.abspath( os.path.join( datadir, fname ) ) )
            group.append( int( match.group( 1 ) ) )
    return zip( *sorted( zip( group, result ) ) ) 
    #return sort( result ), group
    
num_slice, slice_files = get_slice_files( data_dir )

n_slices = len( slice_files )
start_slice = npmin( num_slice )
slice_range = range( start_slice, start_slice + n_slices )
slice_distance = 500 # micrometers


def read_yarn_structure():

    slice_point_list = []
    slice_radius_list = []
    slice_len_list = []

    #cut_off_start = zeros( ( n_slices, ), dtype='int' )
    #cut_off_start[ 1: ] += 0

    for slice_idx, data_file in zip( slice_range, slice_files ):

        #data_file = join( data_dir, 'V1Schnitt%d.txt' % slice_idx )

        print 'reading data_file'

        points = loadtxt( data_file ,
                          skiprows=1,
                          usecols=( 1, 2, 3 ) )

        y = points[ :, 0]
        z = points[ :, 1]
        x = ones_like( y ) * slice_idx * slice_distance
        r = points[ :, 2]

        slice_point_list.append( c_[ x, y, z ] )
        slice_radius_list.append( r )
        slice_len_list.append( points.shape[0] )

    lens_arr = array( slice_len_list )
    print 'slice lens', lens_arr
    offset_arr = cumsum( lens_arr )
    slice_offset_arr = zeros_like( offset_arr )
    slice_offset_arr[1:] = offset_arr[:-1]
    print 'slice offsets', slice_offset_arr

    data_file = join( data_dir, 'connectivity.txt' )
    filam_connect_arr = loadtxt( data_file )
    filam_connect_arr = filam_connect_arr[ npmin( filam_connect_arr, axis=1 ) != -1. ]
    filam_connect_arr = filam_connect_arr[ npmin( filam_connect_arr, axis=1 ) != 0.0 ]
    print filam_connect_arr.shape

    print filam_connect_arr.shape
    #print slice_offset_arr.shape

    fil_map = array( filam_connect_arr + slice_offset_arr, dtype='int' )
    
    points = vstack( slice_point_list )
    radius = hstack( slice_radius_list )

    print points.shape
    print max( fil_map.flatten() )

    p = points[ fil_map.flatten() ]
    r = radius[ fil_map.flatten() ]

    mlab.plot3d( p[:, 0], p[:, 1], p[:, 2], r,
                 tube_radius=20, tube_sides=20, colormap='Spectral' )#
    
    offset = array( [0, 3, 6] )
    cells = array( [10, 4000, 20, 5005, 20, 4080, 4000, 20, 404 ] )

#    line_type = tvtk.Line().cell_type # VTKLine == 10
#    cell_types = array( [line_type] )
#    # Create the array of cells unambiguously.
#    cell_array = tvtk.CellArray()
#    cell_array.set_cells( 3, cells )

    # Now create the UG.
    ug = tvtk.UnstructuredGrid( points=points )
    # Now just set the cell types and reuse the ug locations and cells.
#    ug.set_cells( cell_types, offset, cell_array )
    ug.point_data.scalars = radius
    ug.point_data.scalars.name = 'radius'
    ug.save()
    return ug


# Now view the data.
@mayavi2.standalone
def view( ug ):
    from enthought.mayavi.sources.vtk_data_source import VTKDataSource
    from enthought.mayavi.modules.outline import Outline
    from enthought.mayavi.modules.surface import Surface
    from enthought.mayavi.modules.vectors import Vectors

    mayavi.new_scene()
    src = VTKDataSource( data=ug )
    mayavi.add_source( src )
    s = Surface()
    mayavi.add_module( s )

if __name__ == '__main__':
    ug = read_yarn_structure()
    mlab.show() # view( ug )


