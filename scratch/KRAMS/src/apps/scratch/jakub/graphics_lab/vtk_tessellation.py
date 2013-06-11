'''
Created on Apr 16, 2009

@author: jakub
'''
# Here you can add your own code ...
#!/usr/bin/env python

"""A simple example demonstrating how one can use numpy arrays
transparently with TVTK.

"""
from enthought.traits.api import \
     Array, Bool, Callable, Enum, Float, HasTraits, Interface, implements, \
     Instance, Int, Trait, Str, Enum, Callable, List, TraitDict, Any, \
     on_trait_change, Tuple, WeakRef, Delegate, Property, cached_property

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
#mesh = tvtk.UnstructuredGrid()

class AdaptorCell(tvtk.GenericAdaptorCell):
    
    def __init__(self):
        pass
    
class TesDataSet(tvtk.GenericDataSet):

    points = Array(float)
    
    connectivity = Array(int)
    
    def __init__(self):
        pass
    
class TesCellIterator(tvtk.GenericCellIterator):
    
    def __init__(self):
        pass

mesh = TesDataSet()
mesh.points = points_arr
mesh.connectivity = id_map
iterator = TesCellIterator()
quad = AdaptorCell()




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