
from numpy import array
# tvtk related imports
#
from enthought.tvtk.api import tvtk
from enthought.mayavi.scripts import mayavi2
from enthought.mayavi.sources.vtk_data_source import VTKDataSource
from enthought.mayavi.modules.api import Outline, Surface

@mayavi2.standalone
def view():
    # 'mayavi' is always defined on the interpreter.
    mayavi.new_scene()
    # Make the data and add it to the pipeline.

    height =  1.
    width = 1.

    points =  array([[0,0,0],
                  [width/2.,   0,         0],
                  [width/2.,   height/2., 0],
                  [0,          height/2., 0]] )
    
    faces = array([[0,3,2,1]])

    poly_data = tvtk.PolyData(points = points, polys = faces)
    poly_data.point_data.scalars = array([1.,4.,3.,3])
    poly_data.point_data.scalars.name = 'My spatial function'

    src = VTKDataSource( data = poly_data )
    mayavi.add_source( src )
    # Visualize the data.
    mayavi.add_module(Outline())
    mayavi.add_module(Surface())    
    
view()