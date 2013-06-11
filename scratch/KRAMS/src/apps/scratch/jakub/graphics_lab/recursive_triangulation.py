'''
Created on Jul 1, 2009

@author: jakub
'''
from enthought.tvtk.api import tvtk
# Generate some random points
from numpy import array, sum, ix_, average, int_, vstack

points = tvtk.Points()

points.insert_next_point(0, 0, 0)
points.insert_next_point(0.5, 0, 0)
points.insert_next_point(0.5, 1, 0)
#points.insert_next_point(1, 4, 0)
#points.insert_next_point(0.5, 3.5, 0)
points.insert_next_point(0, 1, 0)
points.insert_next_point(0.25, 0.5, 0)

# Create a polydata with the points we just created.
#profile = tvtk.PolyData(points=points)
profile = tvtk.PolyData(points=points, polys = array([[1,2,3,4]]))

# Perform a 2D Delaunay triangulation on them.
delny = tvtk.Delaunay2D(input=profile, offset = 1.e1)
tri= delny.output
tri.update()
#for i in range(len(tri.polys.data)):
#    print i," ",tri.polys.data.get_value(i)
triangles = array(tri.polys.data,dtype=int_)    
#print "\n"
#for i in range(len(tri.points.data)):
#    print i," ",tri.points.data[i][0]
points = array(tri.points.data)
print "points ", points
#print "triangles ", triangles
ids = (triangles.reshape((triangles.shape[0]/4),4))[:,1:]
print "ids ",ids
#
#weigths = array([[0.5,0.5,0.],[0.,0.5,0.5],[0.5,0.,0.5]])
#
#new_pnts = vstack((average(points[ix_(ids[0])],0,weigths[0]),
#                   average(points[ix_(ids[0])],0,weigths[1]),
#                   average(points[ix_(ids[0])],0,weigths[2])))
#
#print "new points ", new_pnts
#
#new_input = vstack((points[ids[0]],new_pnts))
#
#second = tvtk.PolyData(points=new_input)
#delny = tvtk.Delaunay2D(input=second)
#tri= delny.output
#tri.update()
##for i in range(len(tri.polys.data)):
##    print i," ",tri.polys.data.get_value(i)
#triangles = array(tri.polys.data,dtype=int_)    
##print "\n"
##for i in range(len(tri.points.data)):
##    print i," ",tri.points.data[i][0]
#points = array(tri.points.data)
#print "points ", points
##print "triangles ", triangles
#ids = (triangles.reshape((triangles.shape[0]/4),4))[:,1:]
#print "ids ",ids
