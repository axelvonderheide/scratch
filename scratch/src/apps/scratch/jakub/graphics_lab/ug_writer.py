'''
(C)2004 ORFEUS

This is the default python module file. It has been generated
by the modtree administrator

Import the generated Python module in order to
 - initialize the required modules and
 - load the compiled  module extension library of this module
'''


# Here you can add your own code ...
#!/usr/bin/env python

"""A simple example demonstrating how one can use numpy arrays
transparently with TVTK.

"""

# Author: Prabhu Ramachandran and Eric Jones
# Copyright (c) 2004-2007, Enthought, Inc.
# License: BSD Style.

from enthought.tvtk.api import tvtk
from numpy import array




temperature = tvtk.DoubleArray()
temperature.insert_next_value(0.)
temperature.insert_next_value(0.)
temperature.insert_next_value(0.)
temperature.insert_next_value(0.)

temperature.insert_next_value(20.)
temperature.insert_next_value(0.)
temperature.insert_next_value(0.)
temperature.insert_next_value(0.)

#temp_array.append(temp2)
### TVTK PIPELINE
# create a render window and hand it the renderer
render_window = tvtk.RenderWindow(size=(400,400))

# create a renderer
renderer = tvtk.Renderer(background=(0.839216, 0.839216, 0.839216))
render_window.add_renderer(renderer)

# create interactor and hand it the render window
# This handles mouse interaction with window.
interactor = tvtk.RenderWindowInteractor(render_window=render_window)



points_arr = tvtk.Points()
#points_arr.insert_point(0, -1, -1, 0)
#points_arr.insert_point(1, 1, -1, 0)
#points_arr.insert_point(2, 1, 1, 0)
#points_arr.insert_point(3, -1, 1, 0)
#points_arr.insert_point(0, 0., 0., 0.)
#points_arr.insert_point(1, 0.5, 0., 1.)
#points_arr.insert_point(2, 1., 0., 0.)
#points_arr.insert_point(3, 1., 0.5, 0.)
#points_arr.insert_point(4, 1, 1., 0.)
#points_arr.insert_point(5, 0.5, 1, 1.)
#points_arr.insert_point(6, 0, 1., 0.)
#points_arr.insert_point(7, 0., 0.5, 0.)

points_arr.insert_point(0, 0., 0., 0.)
points_arr.insert_point(1, 1., 0., 0.)
points_arr.insert_point(2, 1., 1., 0.)
points_arr.insert_point(3, 0., 1., 0.)

points_arr.insert_point(4, 0.5, 0., 0.)
points_arr.insert_point(5, 1, 0.5, 0.)
points_arr.insert_point(6, 0.5, 1., 0.)
points_arr.insert_point(7, 0., 0.5, 0.)

#points_arr.insert_point(4, 1.5, 0., 0.)
#points_arr.insert_point(5, 2, 0.5, 0.)
#points_arr.insert_point(6, 1.5, 1., 0.)
#points_arr.insert_point(7, 1., 0.5, 0.)


quad = tvtk.QuadraticQuad()
quad._get_point_ids().set_id(0, 0)
quad._get_point_ids().set_id(1, 1)
quad._get_point_ids().set_id(2, 2)
quad._get_point_ids().set_id(3, 3)
quad._get_point_ids().set_id(4, 4)
quad._get_point_ids().set_id(5, 5)
quad._get_point_ids().set_id(6, 6)
quad._get_point_ids().set_id(7, 7)

id_map = array( [0,1,2,3,4,5,6,7] )
print "type", quad.cell_type

#define array
#polys = tvtk.CellArray()
#connect a cell to the array
#polys.insert_next_cell(quad)
 
# Create a mesh from the data created above.
#mesh = tvtk.PolyData(points = points_arr,polys = polys)
mesh = tvtk.UnstructuredGrid()
#mesh.insert_next_cell(quad.cell_type,quad._get_point_ids() )
mesh.insert_next_cell(23,id_map )
mesh.points = points_arr
mesh.point_data.scalars = temperature

writer = tvtk.UnstructuredGridWriter()
writer.set_input(mesh)
writer.file_name = "test.vtk"
writer.scalars_name = "temperature"
writer.write()

# Set the mapper to scale temperature range
# across the entire range of colors
#mapper = tvtk.PolyDataMapper(input=mesh)
mapper = tvtk.DataSetMapper(input=mesh)


# Create mesh actor for display
actor = tvtk.Actor(mapper=mapper)
actor.property.color=(1, 0, 0)
actor.property.point_size=(200.)
actor.property.line_width=(200.)




# Now add the actors to the renderer and start the interaction.
renderer.add_actor(actor)

interactor.initialize()
interactor.start()
