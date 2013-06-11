'''
Created on Apr 18, 2009

@author: jakub
'''
'''
Created on Mar 6, 2009

@author: jakub
'''

from enthought.tvtk.api import tvtk
from numpy import array,arange, zeros_like, zeros, sqrt, vstack,\
    linspace, diag, dot
import math

a = 2.
K = 1.
E = 100. 

b = 4.
h = 4.
fin_crack = 20
fin_front = 20

def w_func(A):
    cc = 4 * K /math.sqrt(2*math.pi)/E
    return cc * sqrt(a-A[:,0])

def get_w_arr(A):
    B = zeros_like(A)
    #d_idx = A[:,0]<=a
    #B[d_idx,1]=w_func(A[d_idx])
    B[:,1]=w_func(A)
    return B

def s_func(A):
    cc = K /math.sqrt(2*math.pi)
    return cc / sqrt(A[:,0]- a + (b-a)/fin_front)

def get_s_arr(A):
    B = zeros((A.shape[0],1))
    #d_idx = A[:,0]>=a
    #B[d_idx,0]=s_func(A[d_idx])
    B[:,0] = s_func(A)
    return B
    
x_coords_c = linspace(0.,a,fin_crack)#crack points
x_n_coords_c = linspace(a,0.,20)

x_coords_f = linspace(a,b,fin_front)#front points
x_n_coords_f = linspace(b,a, 20)

points_arr_c = zeros((x_coords_c.shape[0],3))#points 3d coords
points_arr_c[:,0] = x_coords_c

points_arr_f = zeros((x_coords_f.shape[0],3))
points_arr_f[:,0] = x_coords_f

corner_pts_c = zeros((x_n_coords_c.shape[0],3))
corner_pts_c[:,0] = x_n_coords_c
corner_pts_c[:,1] = h/2.

corner_pts_f = zeros((x_n_coords_f.shape[0],3))
corner_pts_f[:,0] = x_n_coords_f
corner_pts_f[:,1] = h/2.

displ_c = get_w_arr(points_arr_c)
points_c1 = vstack((points_arr_c,corner_pts_c))
points_f1 = vstack((points_arr_f,corner_pts_f))

displ_c1 = zeros_like(points_c1)
displ_c1[0:fin_crack,:] = displ_c

displ_f1 = zeros_like(points_f1)

stress_c1 = zeros((points_c1.shape[0],1))
stress_f = get_s_arr(points_arr_f)

stress_f1 = zeros((points_f1.shape[0],1))
stress_f1[0:fin_front] = stress_f
    
points1 = vstack((points_c1,points_f1))
displ_1 = vstack((displ_c1,displ_f1))
stress_1 = vstack((stress_c1,stress_f1))
id_map1 = arange( points_c1.shape[0] )
id_map2 = arange( points_f1.shape[0] ) + id_map1.shape[0]
id_map3 = id_map1 + id_map1.shape[0] + id_map2.shape[0]
id_map4 = id_map2 + id_map1.shape[0] + id_map2.shape[0]

t_mtx = diag([1.,-1.,1.]) 

displ_2 = dot(displ_1,t_mtx)
displ_all = vstack((displ_1,displ_2))

points2 = dot(points1,t_mtx)
points_all = vstack((points1,points2))

stress_all = vstack((stress_1,stress_1))
stress = tvtk.DoubleArray(name = 'stress')
stress.from_array(stress_all)

displ = tvtk.DoubleArray(name = 'displ')
displ.from_array(displ_all)


mesh = tvtk.UnstructuredGrid()
mesh.points = points_all
#mesh.insert_next_cell(quad.cell_type,quad._get_point_ids() )
mesh.insert_next_cell(7,id_map1 )
mesh.insert_next_cell(7,id_map2 )
#mesh.insert_next_cell(7,id_map3 )
#mesh.insert_next_cell(7,id_map4 )


#mesh.point_data.scalars = temperature

mesh.point_data.add_array(stress)
mesh.point_data.add_array(displ)


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
    from enthought.mayavi.filters.api import WarpVector, ExtractTensorComponents,\
        WarpScalar, ExtractEdges, Delaunay2D

    mayavi.new_scene()
    src = VTKDataSource(data = mesh)
    mayavi.add_source(src) 
    mayavi.add_filter(ExtractEdges())
    mayavi.add_filter(Delaunay2D())
    mayavi.add_filter( WarpVector())
    mayavi.add_filter( WarpScalar())
    mayavi.add_module(Surface())


if __name__ == '__main__':
    view()
