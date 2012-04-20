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
temperature.insert_next_value(20.)
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

_node_coord_map = array( [[-1., -1.],
                          [ 1., -1.],
                          [ 1.,  1.],
                          [-1.,  1.]] )


def get_N_geo_mtx(r_pnt):
    '''
    Return the value of shape functions for the specified local coordinate r
    '''
    cx = _node_coord_map
    Nr = array( [[1/4.*(1 + r_pnt[0]*cx[i,0])*(1 + r_pnt[1]*cx[i,1]) 
                  for i in range(0,4) ]] )
    return Nr

def get_dNr_geo_mtx( r_pnt):
    '''
    Return the matrix of shape function derivatives.
    Used for the conrcution of the Jacobi matrix.

    @TODO - the B matrix is used
    just for uniaxial bar here with a trivial differential
    operator.
    '''
    cx = _node_coord_map
    dNr_geo = array( [[ 1/4.*cx[i,0]*(1 + r_pnt[1]*cx[i,1]) for i in range(0,4) ],
                      [ 1/4.*cx[i,1]*(1 + r_pnt[0]*cx[i,0]) for i in range(0,4) ]])        
    return dNr_geo



#    def interpolation_derivs(self, *args):
#        """
#        V.interpolation_derivs((float, float, float), (float, float, float, float, float, float, float, float))
#        
#         Pixel specific methods.
#        """
#        ret = self._wrap_call(self._vtk_obj.InterpolationDerivs, *args)
#        return ret
#
#    def interpolation_functions(self, *args):
#        """
#        V.interpolation_functions((float, float, float), (float, float, float, float))
#        
#         Pixel specific methods.
#        """
#        ret = self._wrap_call(self._vtk_obj.InterpolationFunctions, *args)
#        return ret


points_arr = tvtk.Points()
#points_arr.insert_point(0, -1, -1, 0)
#points_arr.insert_point(1, 1, -1, 0)
#points_arr.insert_point(2, 1, 1, 0)
#points_arr.insert_point(3, -1, 1, 0)
points_arr.insert_point(0, 0, 0, 0)
points_arr.insert_point(1, 1, 0, 0)
points_arr.insert_point(2, 1, 1, 0)
points_arr.insert_point(3, 0, 1, 0)



quad = tvtk.Pixel()
quad._get_point_ids().set_id(0, 0)
quad._get_point_ids().set_id(1, 1)
quad._get_point_ids().set_id(2, 2)
quad._get_point_ids().set_id(3, 3)


#aaa = quad.interpolation_functions((0.1,0.1,0.),(0.,0.,1.,0.))
#print 'aaa', aaa


#define array
polys = tvtk.CellArray()

#connect a cell to the array
polys.insert_next_cell(quad)
 






# Create a mesh from the data created above.
mesh = tvtk.PolyData(points = points_arr,polys = polys)
mesh.point_data.scalars = temperature
# Set the mapper to scale temperature range
# across the entire range of colors
mapper = tvtk.PolyDataMapper(input=mesh)


# Create mesh actor for display
actor = tvtk.Actor(mapper=mapper)
actor.property.color=(1, 0, 0)
actor.property.point_size=(200.)
actor.property.line_width=(200.)


# Now add the actors to the renderer and start the interaction.
renderer.add_actor(actor)

interactor.initialize()
interactor.start()
