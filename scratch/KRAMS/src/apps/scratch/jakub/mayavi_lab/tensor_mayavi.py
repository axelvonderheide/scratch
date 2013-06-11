'''
Created on Mar 6, 2009

@author: jakub
'''

from enthought.tvtk.api import tvtk
from numpy import array

temperature = tvtk.DoubleArray()
temperature.insert_next_value(0.)
temperature.insert_next_value(0.)
temperature.insert_next_value(2.)
temperature.insert_next_value(2.)

strain = array([[0.,0,0,0,2,0,0,0,1],
                [2,0,0,0,1,0,0,0,0],
                [1,0,0,0,1,0,0,0,2],
                [1,0,0,0,1,0,0,0,0]])

stress = array([[0.,0,0,0,2,0,0,0,1],
                [8.,0,0,0,1,0,0,0,0],
                [9.,0,0,0,1,0,0,0,2],
                [6.,0,0,0,1,0,0,0,0]])

displ1 = array([[0.,0,0],
                [8.,0,0],
                [9.,0,0],
                [6.,0,0]])

displ2 = array([[1.,0,0],
                [2.,0,0],
                [3.,0,0],
                [4.,0,0]])

temperature = array([[1.],
                [0.],
                [0.],
                [0.]])

str1 = tvtk.DoubleArray(name = 'first')
str1.from_array(strain)


str2 = tvtk.DoubleArray(name = 'second')
str2.from_array(stress)

temp = tvtk.DoubleArray(name = 'temp')
temp.from_array(temperature)

d1 = tvtk.DoubleArray(name = 'd1')
d1.from_array(displ1)

d2 = tvtk.DoubleArray(name = 'd2')
d2.from_array(displ2)

points_arr = [[0.,0.,0.],
              [1.,0.,0.],
              [0.,1.,0.],
              [1.,1.,0.]
              ]
id_map = array( [0,1,2,3] )
polys = array( [[0,1,2,3]] )




mesh = tvtk.UnstructuredGrid()
#mesh.insert_next_cell(quad.cell_type,quad._get_point_ids() )
mesh.insert_next_cell(8,id_map )
mesh.points = points_arr
#mesh.point_data.scalars = temperature

mesh.point_data.add_array(str1)
mesh.point_data.add_array(str2)
#mesh.point_data.tensors = str1
#mesh.point_data.tensors = str2
mesh.point_data.add_array(temp)
mesh.point_data.add_array(d1)
mesh.point_data.add_array(d2)

#mesh.point_data.set_active_tensors('second')
#print "str1", mesh.point_data.tensors
#test

#print "str2", str2
#print mesh.is_array_an_attribute

#mesh.point_data.tensors = mesh.field_data.get_array('first')
#mesh.point_data.tensors = str1

#mesh.point_data.set_active_attribute(1,4)
#print 'test',mesh.point_data.is_array_an_attribute(1)

from enthought.mayavi.scripts import mayavi2
# Now view the data.
@mayavi2.standalone
def view():
    from enthought.mayavi.sources.vtk_data_source import VTKDataSource
    from enthought.mayavi.modules.outline import Outline
    from enthought.mayavi.modules.surface import Surface
    from enthought.mayavi.modules.vectors import Vectors
    from enthought.mayavi.modules.api import TensorGlyph
    from enthought.mayavi.filters.api import WarpVector, ExtractTensorComponents

    mayavi.new_scene()
    
    wv = WarpVector()
    tc = ExtractTensorComponents()
    tc2 = ExtractTensorComponents()
    # The single type one
    src = VTKDataSource(data = mesh)
    print "str1", mesh.point_data.is_array_an_attribute(0)
    print "id ", mesh.point_data.has_array("d2")
    mayavi.add_source(src) 
    mayavi.add_filter(wv)
    src.add_child(tc)
    tc.add_module(Surface())
    #src.add_child(Surface())
    
    mesh.point_data.set_active_tensors('first')
    src2 = VTKDataSource(data = mesh)
    
    mayavi.add_source(src2) 
    #tg = TensorGlyph()
    #tg.glyph.glyph.scale_factor = 100.
    #mayavi.add_module(tg)
    #mayavi.add_module(Surface())
    #mayavi.add_module(Vectors())
    src2.add_child(tc2)
    tc2.add_module(Surface())
    tc2.scene.z_plus_view()

if __name__ == '__main__':
    view()
