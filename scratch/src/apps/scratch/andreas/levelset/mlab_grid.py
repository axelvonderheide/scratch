'''
Created on Sep 2, 2009

@author: jakub
'''
import numpy
from enthought.mayavi.mlab import surf, colorbar, show, mesh

@show
def gsurf(): 
    """Test contour_surf on regularly spaced co-ordinates like MayaVi."""
    def f(x, y):
        return x**2+y**2 - 25.

    x, y = numpy.mgrid[-10.:10.:0.25, -10.:10.0:0.25]
    s = surf(x, y, f, warp_scale = 0.)
    colorbar(s)
    return s



if __name__ == '__main__':
    gsurf()
