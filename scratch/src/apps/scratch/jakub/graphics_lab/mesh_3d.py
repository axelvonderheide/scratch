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

### DATA
points = array([[0,0,0],
                [2,0,0],
                [2,2,0],
                [0,2,0],
                [0,0,2],
                [2,0,2],
                [2,2,2],
                [0,2,2]], 'f')

#vertices = array([[0],
#                  [1],
#                  [2],
#                  [3],
#                  [4]])

#lines = array([[0,1],
#               [1,2],
#               [2,3],
#               [3,0]])

#triangles = array([[0,1,2],
#                   [0,2,3]])


triangles = tvtk.CellArray()

hex= tvtk.Hexahedron()

triangles.insert_next_cell(hex)
#triangles.insert_cell_point(0)
#triangles.insert_cell_point(1)
#triangles.insert_cell_point(2)
#triangles.insert_cell_point(3)
#triangles.insert_cell_point(4)
#triangles.insert_cell_point(5)
#triangles.insert_cell_point(6)
#triangles.insert_cell_point(7)

#temperature = [[0.,20.,20.,0.]]
#temperature = tvtk.DoubleArray()
#temperature.insert_next_value(0.)
#temperature.insert_next_value(20.)
#temperature.insert_next_value(20.)
#temperature.insert_next_value(0.)
#temperature.insert_next_value(0.)
#temperature.insert_next_value(0.)
#temperature.insert_next_value(20.)
#temperature.insert_next_value(0.)
#
#temp_array=[]
#temp1 = tvtk.DoubleArray()
#temp1.insert_next_value(0.)
#temp1.insert_next_value(20.)
#temp1.insert_next_value(20.)
#temp1.insert_next_value(0.)
#temp_array.append(temp1)

#temp2 = tvtk.DoubleArray()
#temp2.insert_next_value(0.)
#temp2.insert_next_value(20.)
#temp2.insert_next_value(20.)
#temp2.insert_next_value(0.)
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




#points_arr = tvtk.Points()
#points_arr.insert_next_point(0, 0, 0)
#points_arr.insert_next_point(1, 0, 0)
#points_arr.insert_next_point(1, 1, 0)
#points_arr.insert_next_point(0, 1, 0)

## ids = tvtk.IdList()
## ids.insert_next_id(0)
## ids.insert_next_id(1)
## ids.insert_next_id(2)
## ids.insert_next_id(3)



#polys = tvtk.CellArray()
#polys.insert_next_cell(3)
#polys.insert_cell_point(0)
#polys.insert_cell_point(1)
#polys.insert_cell_point(2)
#polys.insert_next_cell(3)
#polys.insert_cell_point(0)
#polys.insert_cell_point(2)
#polys.insert_cell_point(3)


# Create a mesh from the data created above.
mesh = tvtk.PolyData(points = points,polys = triangles)
#v_mesh = tvtk.PolyData(points = points, verts = vertices )
#c_mesh = tvtk.PolyData(points = points,lines = lines)
#p_mesh = tvtk.PolyData(points = points,polys = triangles)
#p_mesh.point_data.scalars = temp_array


# Set the mapper to scale temperature range
# across the entire range of colors
mapper = tvtk.PolyDataMapper(input=mesh)
#mapper.interpolate_scalars_before_mapping=False
#v_mapper = tvtk.PolyDataMapper(input=v_mesh)
#c_mapper = tvtk.PolyDataMapper(input=c_mesh)
#p_mapper = tvtk.PolyDataMapper(input=p_mesh)
#p_mapper.interpolate_scalars_before_mapping=False


# Create mesh actor for display
actor = tvtk.Actor(mapper=mapper)
actor.property.color=(1, 0, 0)
actor.property.point_size=(200.)
actor.property.line_width=(200.)


# Now add the actors to the renderer and start the interaction.
renderer.add_actor(actor)
#renderer.add_actor(v_actor)
interactor.initialize()
interactor.start()
