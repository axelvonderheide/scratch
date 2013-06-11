'''
Created on Jul 1, 2009

@author: jakub
'''
from numpy import array
from matplotlib.delaunay.triangulate import Triangulation

x = array([0.5,0.,1.,1.,0])
y = array([0.5,0.,0.,1.,1.])

tri = Triangulation(x,y)
print tri.circumcenters 

print tri.edge_db 

print tri.triangle_nodes

print tri.triangle_neighbors