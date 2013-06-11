'''
Created on Dec 22, 2009

@author: alexander
'''

from numpy import array, loadtxt, vstack, c_

def read_rsurface( filename ):
    '''Read the rhino data format.
    Extract solely the vertex positions.
    '''
    geo_file = file( filename )
    vertices_list = []
    for line in geo_file:
        if line[0] == "v" and line[0:2] != "vp":
            coords = line.split()
            x, y, z = map( float, coords[1:] )
            vertices_list.append( [x,y,z] )
    
    vertices_arr = array(vertices_list)
    return vertices_arr

def normalize_rsurface( v_arr ):
    vx_arr, vy_arr, vz_arr = v_arr[:,0], v_arr[:,1], v_arr[:,2]
    xmin, xmax = min( vx_arr ), max( vx_arr )
    ymin, ymax = min( vy_arr ), max( vy_arr )
    zmin, zmax = min( vz_arr ), max( vz_arr )

    xrange = xmax - xmin 
    yrange = ymax - ymin
    zrange = zmax - zmin

    vx_arr = ( vx_arr - xmin ) * 1 / xrange
    vy_arr = ( vy_arr - ymin ) * 1 / yrange
    vz_arr = ( vz_arr - zmin ) * 1 / zrange

    return c_[ vx_arr, vy_arr, vz_arr ]
    
if __name__ == '__main__':
    #geo_file = file( 'lowerface.robj' )
    geo_filename = 'upperface.robj'

    v_arr = read_rsurface( geo_filename )
    print v_arr
    print v_arr.shape
    print '\n'
    
    v_arr = normalize_rsurface( v_arr )
    print 'normalized'
    print v_arr


    
