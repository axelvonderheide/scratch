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
points = array([[2,0,0],
                [4,0,0],
                [4,2,0],
                [2,2,0]], 'f')

vertices = array([[0],
                  [1],
                  [2],
                  [3]])

lines = array([[0,1],
               [1,2],
               [2,3],
               [3,0]])

#triangles = array([[0,1,2],
#                   [0,2,3]])

triangles = array([[0,1,2,3]])
#triangles = tvtk.CellArray()
#triangles.insert_next_cell(4)
#triangles.insert_cell_point(0)
#triangles.insert_cell_point(1)
#triangles.insert_cell_point(2)
#triangles.insert_cell_point(3)

#temperature = [[0.,20.,20.,0.]]
temperature = tvtk.DoubleArray()
temperature.insert_next_value(0.)
temperature.insert_next_value(20.)
temperature.insert_next_value(20.)
temperature.insert_next_value(0.)




deformation = tvtk.DoubleArray()
deformation.number_of_components = 3
deformation.number_of_tuples = 4
deformation.set_tuple3(0,0.,0.,0)
deformation.set_tuple3(1,20.,-5.,0.)
deformation.set_tuple3(2,15.,3.,0.)
deformation.set_tuple3(3,-4.,2.,0)
#temp_array.append(deformation)

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
points_arr.insert_next_point(0, 0, 0)
points_arr.insert_next_point(1, 0, 0)
points_arr.insert_next_point(1, 1, 0)
points_arr.insert_next_point(0, 1, 0)

## ids = tvtk.IdList()
## ids.insert_next_id(0)
## ids.insert_next_id(1)
## ids.insert_next_id(2)
## ids.insert_next_id(3)



polys = tvtk.CellArray()
polys.insert_next_cell(3)
polys.insert_cell_point(0)
polys.insert_cell_point(1)
polys.insert_cell_point(2)
polys.insert_next_cell(3)
polys.insert_cell_point(0)
polys.insert_cell_point(2)
polys.insert_cell_point(3)


# Create a mesh from the data created above.
mesh = tvtk.PolyData(points = points,polys = polys)
v_mesh = tvtk.PolyData(points = points, verts = vertices )
v_mesh.point_data.vectors = deformation
c_mesh = tvtk.PolyData(points = points,lines = lines)
c_mesh.point_data.vectors = deformation
p_mesh = tvtk.PolyData(points = points,lines = lines,polys = triangles)
p_mesh.point_data.scalars = temperature
p_mesh.point_data.vectors = deformation

warp_v = tvtk.WarpVector(input=v_mesh, scale_factor = .2)
warp_c = tvtk.WarpVector(input=c_mesh, scale_factor = .0)
warp_p = tvtk.WarpVector(input=p_mesh, scale_factor = .2)

# Set the mapper to scale temperature range
# across the entire range of colors
mapper = tvtk.PolyDataMapper(input=mesh)
mapper.interpolate_scalars_before_mapping=False
v_mapper = tvtk.PolyDataMapper(input=warp_v.output)
#v_mapper = tvtk.PolyDataMapper(input=v_mesh)
#c_mapper = tvtk.PolyDataMapper(input=c_mesh)
c_mapper = tvtk.PolyDataMapper(input=warp_c.output)
#p_mapper = tvtk.PolyDataMapper(input=p_mesh)
p_mapper = tvtk.PolyDataMapper(input=warp_p.output)


sphere = tvtk.SphereSource(radius = 0.05,phi_resolution = 20, theta_resolution=20)

ball = tvtk.Glyph3D(input = v_mesh,source = sphere.output)
#vertices.SetInputConnection(ThresholdIn.GetOutputPort())
#vertices.SetSource(Sphere.GetOutput())

ball_mapper = tvtk.PolyDataMapper(input= ball.output)

ball_actor = tvtk.Actor(mapper=ball_mapper)
ball_actor.property.color=(0, 0, 0)
ball_actor.property.point_size=(200.)

renderer.add_actor(ball_actor)

# Create mesh actor for display
actor = tvtk.Actor(mapper=mapper)
actor.property.color=(1, 0, 0)
actor.property.point_size=(200.)
actor.property.line_width=(200.)

v_actor = tvtk.Actor(mapper=v_mapper)
v_actor.property.color=(0, 0, 0)
v_actor.property.point_size=(5.)

c_actor = tvtk.Actor(mapper=c_mapper)
c_actor.property.color=(0.27451, 0.509804, 0.705882)
c_actor.property.line_width=(4.)

p_actor = tvtk.Actor(mapper=p_mapper)
p_actor.property.color=(2./3, 2./3, 2./3)



text_actor =tvtk.TextActor(position=(200,100))
text_actor.scaled_text= True
text_actor.input = "here"
#text_actor.display_position(2,4)
renderer.add_actor2d(text_actor)

atext = tvtk.VectorText()
atext.text="Origin"
textMapper = tvtk.PolyDataMapper()
textMapper.input_connection=atext.output_port

textActor = tvtk.Follower(position=(0, -0.1, 0.1), mapper=textMapper)
textActor.scale=(0.2, 0.2, 0.2)
textActor.property.color=(0, 0, 0)
#textActor.add_position=(0, -0.1, 0)
renderer.add_actor(textActor)

# Now add the actors to the renderer and start the interaction.
renderer.add_actor(actor)

renderer.add_actor(c_actor)
renderer.add_actor(p_actor)
#renderer.add_actor(v_actor)
interactor.initialize()
interactor.start()
