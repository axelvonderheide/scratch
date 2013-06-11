'''
Created on 08.05.2009

@author: alexander
'''

from numpy import *

shape = (40,10)
n_nodes = (shape[0]+1, shape[1]+1)
print n_nodes
m = mgrid[:n_nodes[0],:n_nodes[1]]
#print m.shape

points = m.flatten()

print points.shape
n_nodes_xy = n_nodes[0]*n_nodes[1]
x = points[:n_nodes_xy]
y = points[n_nodes_xy:]
points = c_[x,y]
print 'points: \n', points
print 'n_nodes[0]: \n', n_nodes[0]
print 'n_nodes[1]: \n', n_nodes[1]

def arch_2d( points ):
    x = points[:,0]
    print 'x: \n', x
    y = points[:,1]
    print 'y: \n', y

    R_a = 1.395
    beta = 100.18 * pi/180
    print 'beta: \n', beta

    D_top = 0.03
    D_bottom = 0.04
    D_delta = D_bottom - D_top

    no_x = arange(n_nodes_xy)/n_nodes[1]
    print 'no_x: \n', no_x
    D_var_left = D_bottom - no_x/(n_nodes[0]-1.) * D_delta * 2
    print 'D_var_left: \n', D_var_left
    D_var_right = D_bottom - ((n_nodes[0]-1.) - no_x)/(n_nodes[0]-1) * D_delta * 2
    print 'D_var_right: \n', D_var_right
    D_var = hstack((D_var_left[:n_nodes_xy/2],D_var_right[n_nodes_xy/2:]))
    print 'D_var: \n', D_var
    
    R_i = R_a - D_bottom
    print 'R_i: \n', R_i
    alpha = (pi - beta )/2
    print 'alpha: \n', alpha
    
    phi = alpha + x / (x[-1]-x[0]) * beta
    print 'phi: \n', phi.shape[0]

    r = R_i + y / (y[-1]-y[0]) * D_var
    x,y = - r * cos( phi ), r * sin( phi ) 
    return c_[ x,y ]

points = arch_2d( points )




#shape = (40,10)
#R_a = 1.395
#beta = 100.18 * pi/180
#D_top = 0.03
#D_bottom = 0.04
#
#D_delta = D_bottom - D_top
#n_nodes = (shape[0]+1, shape[1]+1)
#n_nodes_xy = n_nodes[0]*n_nodes[1]
#no_x = arange(n_nodes_xy)/n_nodes[1]
#D_var_left = D_bottom - no_x/(n_nodes[0]-1.) * D_delta * 2
#D_var_right = D_bottom - ((n_nodes[0]-1.) - no_x)/(n_nodes[0]-1) * D_delta * 2
#D_var = hstack((D_var_left[:n_nodes_xy/2],D_var_right[n_nodes_xy/2:]))
#R_i = R_a - D_bottom
#alpha = (pi - beta )/2
#phi = alpha + x / (x[-1]-x[0]) * beta
#r = R_i + y / (y[-1]-y[0]) * D_var
#x,y = - r * cos( phi ), r * sin( phi ) 


