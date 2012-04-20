"""A MayaVi example of how to generate a polydata dataset using numpy
arrays.  The example is similar to tvtk/examples/tiny_mesh.py.  Also
shown is a way to visualize this data with mayavi2.  The script can be
run like so:

  $ mayavi2 -x polydata.py

It can be alternatively run as:

  $ python polydata.py
  
Author: Prabhu Ramachandran <prabhu at aero dot iitb dot ac dot in>

Copyright (c) 2007, Enthought, Inc.
License: BSD style.
"""

# The following *optional* two lines allow a user to call this script
# as either `python script.py` or `mayavi2 script.py`.  These two
# lines must be placed before any other mayavi imports.
from enthought.mayavi.scripts import mayavi2
mayavi2.standalone(globals())

from numpy import array
from enthought.tvtk.api import tvtk

# The numpy array data.
points = array([[0,0,0], [1,0,0], [0,1,0], [0,0,1]], 'f')
triangles = array([[0,1,3], [0,3,2], [1,2,3], [0,2,1]])
temperature = array([10., 20., 30., 40.])

# The TVTK dataset.
mesh = tvtk.PolyData(points=points, polys=triangles)
mesh.point_data.scalars = temperature
mesh.point_data.scalars.name = 'Temperature'

# Uncomment the next two lines to save the dataset to a VTK XML file.
#w = tvtk.XMLPolyDataWriter(input=mesh, file_name='polydata.vtp')
#w.write()

# Now view the data.
def view():
    from enthought.mayavi.sources.vtk_data_source import VTKDataSource
    from enthought.mayavi.modules.surface import Surface

    mayavi2.new_scene()
    src = VTKDataSource(data = mesh)
    mayavi2.add_source(src)
    s = Surface()
    mayavi2.add_module(s)

view()
