'''
Created on Mar 6, 2009

@author: jakub
'''
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

from enthought.tvtk.api import tvtk
from numpy import array

temperature = tvtk.DoubleArray()
temperature.insert_next_value(0.)
temperature.insert_next_value(0.)
temperature.insert_next_value(2.)
temperature.insert_next_value(2.)

#strain = tvtk.DoubleArray()
#strain.insert_next_value(0.)
#strain.insert_next_value(0.)
#strain.insert_next_value(2.)
#strain.insert_next_value(2.)

strain = array([[1,0,0,0,1,0,0,0,1],
                [1,0,0,0,2,0,0,0,1],
                [1,0,0,0,4,0,0,0,1],
                [1,0,0,0,7,0,0,0,1]])

#
#temperature.insert_next_value(0.)
#temperature.insert_next_value(0.)
#temperature.insert_next_value(0.)
#temperature.insert_next_value(0.)
#
#temperature.insert_next_value(0.)
#temperature.insert_next_value(0.)
#temperature.insert_next_value(0.)
#temperature.insert_next_value(0.)
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
              [1.,1.,0.],
              [0.,1.,0.]
              ]
id_map = array( [0,1,2,3] )
polys = array( [[0,1,2,3]] )


#temperature= array( [1.,2.,3.,4.] )
# Create a mesh from the data created above.

mesh = tvtk.UnstructuredGrid()
#mesh.insert_next_cell(quad.cell_type,quad._get_point_ids() )
mesh.insert_next_cell(7,id_map )
mesh.points = points_arr
mesh.point_data.scalars = temperature
mesh.point_data.tensors = strain
# Set the mapper to scale temperature range
# across the entire range of colors
#mapper = tvtk.PolyDataMapper(input=mesh)

#mesh_poly = tvtk.PolyData(points = points_arr,polys = polys)
#mesh_poly.points = points_arr
#mesh_poly.point_data.scalars = temperature
#mesh_poly.point_data.tensors = strain

mesh_poly = tvtk.UnstructuredGridToPolyDataFilter(input = mesh)

tg = tvtk.TensorGlyph(source = mesh_poly,input = mesh_poly)

mapper = tvtk.DataSetMapper(input=tg.output)



# Create mesh actor for display
actor = tvtk.Actor(mapper=mapper)
#actor.property.color=(1, 0, 0)
#actor.property.point_size=(200.)
#actor.property.line_width=(200.)



# Now add the actors to the renderer and start the interaction.
renderer.add_actor(actor)

interactor.initialize()
interactor.start()
