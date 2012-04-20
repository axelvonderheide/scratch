'''
transforms spherical grid (phi, theta, r) to orthogonal coordinates (X,Y,Z) 
'''
from numpy import array, concatenate, arange, transpose, sin, cos,\
                linspace, pi, mgrid, ravel, tile, ones, column_stack

def spheric(r,pt):
    phi = 360 # rotation in XY plane, X axis = 0
    theta = 360 # rotation in XZ plane, Z axis = 0
    grid = mgrid[0:phi/180*pi:phi/180*pi/pt,0:phi/180*pi:phi/180*pi/pt]
    A = column_stack((tile(grid[0,:,0],pt),ravel(grid[0]),ones(pt**2)*r))    
    XYZ = transpose(array([A[:,2]*cos(A[:,0])*sin(A[:,1]),A[:,2]*sin(A[:,0])*sin(A[:,1]),A[:,2]*cos(A[:,1])]))
    return XYZ
     
         