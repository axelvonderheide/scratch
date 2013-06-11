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
temperature.insert_next_value(20.)
temperature.insert_next_value(0.)
temperature.insert_next_value(0.)
temperature.insert_next_value(0.)

temperature.insert_next_value(0.)
temperature.insert_next_value(0.)
temperature.insert_next_value(0.)
temperature.insert_next_value(0.)

temperature.insert_next_value(0.)
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


points_arr = [[0.,0.,0.],
              [1.,0.,0.],
              [2.,0.,0.],
              [3.,0.,0.],
              [3.,1.,0.],
              [3.,2.,0.],
              [3.,3.,0.],
              [2.,3.,0.],
              [ 1.,3.,0.],
              [ 0.,3.,0.],
              [0.,2.,0.],
              [0.,1.,0.]
              ]
id_map = array( [0,1,2,3,4,5,6,7,8,9,10,11] )



#temperature= array( [0,1,2,3,4,5,6,7,8,9,10,11] )
# Create a mesh from the data created above.
#mesh = tvtk.PolyData(points = points_arr,polys = polys)
mesh = tvtk.UnstructuredGrid()
#mesh.insert_next_cell(quad.cell_type,quad._get_point_ids() )
mesh.insert_next_cell(7,id_map )
mesh.points = points_arr
mesh.point_data.scalars = temperature
# Set the mapper to scale temperature range
# across the entire range of colors
#mapper = tvtk.PolyDataMapper(input=mesh)
mapper = tvtk.DataSetMapper(input=mesh)


# Create mesh actor for display
actor = tvtk.Actor(mapper=mapper)
#actor.property.color=(1, 0, 0)
actor.property.point_size=(200.)
actor.property.line_width=(200.)

# Now add the actors to the renderer and start the interaction.
renderer.add_actor(actor)

interactor.initialize()
interactor.start()
