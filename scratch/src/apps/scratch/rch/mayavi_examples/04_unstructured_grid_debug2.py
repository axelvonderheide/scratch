"""A MayaVi example of how to generate an unstructured grid dataset
using numpy arrays.  Also shown is a way to visualize this data with
mayavi2.  The script can be run like so:

  $ mayavi2 -x unstructured_grid.py

Alternatively, it can be run as:

  $ python unstructured_grid.py
  
Author: Prabhu Ramachandran <prabhu at aero dot iitb dot ac dot in>

Copyright (c) 2007, Enthought, Inc.
License: BSD style.
"""

from numpy import array, arange, random
from enthought.tvtk.api import tvtk
from enthought.mayavi.scripts import mayavi2

    
def mixed_type_ug():
    """A slightly more complex example of how to generate an
    unstructured grid with different cell types.  Returns a created
    unstructured grid.
    """
    points = array( [[ 0.,  0.,  0.],
                     [ 2.,  0.,  0.,],
                     [ 2.,  2.,  0.,],
                     [ 0.,  2.,  0.,],
                     [ 0.,  4.,  0.,],
                     [ 2.,  4.,  0.,],
                     [ 2.,  6.,  0.,],
                     [ 0.,  6.,  0.,],
                     [ 0.,  2.,  0.,],
                     [ 2.,  2.,  0.,],
                     [ 2.,  4.,  0.,],
                     [ 0.,  4.,  0.,]], dtype = 'float_' )    
    # The cells
    cells = array( [4,  0,  1,  2,  3,  
                    4,  4,  5,  6,  7,  
                    4,  8,  9, 10, 11], dtype = 'int_')
    # The offsets for the cells, i.e. the indices where the cells
    # start.
    offset = array([0, 5, 10])
    cell_types = array([9, 9, 9])
    # Create the array of cells unambiguously.
    cell_array = tvtk.CellArray()
    cell_array.set_cells(3, cells)
    # Now create the UG.
    ug = tvtk.UnstructuredGrid(points=points)
    # Now just set the cell types and reuse the ug locations and cells.
    ug.set_cells(cell_types, offset, cell_array)
    return ug

    
# ----------------------------------------------------------------------
# Create the unstructured grids and assign scalars and vectors.

ug2 = mixed_type_ug()
temperature = array( [ 2.14285714,
                        -2.14285714,
                        -2.14285714,
                        2.14285714, 
                        0.42857143,
                        -0.42857143,
                        -0.42857143, 
                        0.42857143,
                        1.28571429, 
                        -1.28571429,
                         -1.28571429,
                           1.28571429
                            ], dtype = 'float_' )
ug2.point_data.scalars = temperature
ug2.point_data.scalars.name = 'temperature'

# Uncomment this to save the file to a VTK XML file.
#save_xml(ug2, 'file.vtu')

# Now view the data.
@mayavi2.standalone
def view():
    from enthought.mayavi.sources.vtk_data_source import VTKDataSource
    from enthought.mayavi.modules.outline import Outline
    from enthought.mayavi.modules.surface import Surface
    from enthought.mayavi.modules.vectors import Vectors

    mayavi.new_scene()

    # Mixed types.
    src = VTKDataSource(data = ug2)
    mayavi.add_source(src) 
    mayavi.add_module(Outline())
    mayavi.add_module(Surface())

if __name__ == '__main__':
    view()
