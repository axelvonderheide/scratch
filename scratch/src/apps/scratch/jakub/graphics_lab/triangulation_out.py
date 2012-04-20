from enthought.tvtk.api import tvtk
# Generate some random points
from numpy import array, sum, ix_, average, int_

points = tvtk.Points()
points.insert_next_point(-1, -1, 0)
points.insert_next_point(1, -1, 0)
points.insert_next_point(1, 0, 0)
#points.insert_next_point(1, 4, 0)
#points.insert_next_point(0.5, 3.5, 0)
points.insert_next_point(-1, 0, 0)

# Create a polydata with the points we just created.
profile = tvtk.PolyData(points=points)


# Perform a 2D Delaunay triangulation on them.
delny = tvtk.Delaunay2D(input=profile)
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

def get_gp_list(int_order):
    gps=[]
    if int_order == 1:
        for id in ids:
            gp=[average(points[ix_(id)],0),1.]
            #print "gp ",gp
            gps.append(gp)
    elif int_order == 2:    
        raise NotImplementedError
    elif int_order == 3:    
        weigths = array([[0.6,0.2,0.2],[0.2,0.6,0.2],[0.2,0.2,0.6]])
        for id in ids:
            gp =[average(points[ix_(id)],0),-0.5625],\
                [average(points[ix_(id)],0, weigths[0]),0.52083333333333337],\
                [average(points[ix_(id)],0, weigths[1]),0.52083333333333337],\
                [average(points[ix_(id)],0, weigths[2]),0.52083333333333337]

            gps.append(gp)
    elif int_order == 4:    
        raise NotImplementedError
    elif int_order == 5:    
        weigths = array([[0.0597158717, 0.4701420641, 0.4701420641],\
                         [0.4701420641, 0.0597158717, 0.4701420641],\
                         [0.4701420641, 0.4701420641, 0.0597158717],\
                         [0.7974269853, 0.1012865073, 0.1012865073],\
                         [0.1012865073, 0.7974269853, 0.1012865073],\
                         [0.1012865073, 0.1012865073, 0.7974269853]])
        for id in ids:
            weigts_sum = False#for debug
            gp =[average(points[ix_(id)],0),0.225],\
                [average(points[ix_(id)],0, weigths[0],weigts_sum),0.1323941527],\
                [average(points[ix_(id)],0, weigths[1],weigts_sum),0.1323941527],\
                [average(points[ix_(id)],0, weigths[2],weigts_sum),0.1323941527],\
                [average(points[ix_(id)],0, weigths[3],weigts_sum),0.1259391805],\
                [average(points[ix_(id)],0, weigths[4],weigts_sum),0.1259391805],\
                [average(points[ix_(id)],0, weigths[5],weigts_sum),0.1259391805]

            gps.append(gp)
    else:    
        raise NotImplementedError
    
    return gps
        
        

 #----------------------- example --------------------
if __name__ == '__main__':
    gpl = get_gp_list(5)
    print "gpl ",gpl       