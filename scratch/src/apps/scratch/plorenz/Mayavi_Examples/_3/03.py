import numpy
from enthought.mayavi.mlab import *

def test_points3d():
    t = numpy.linspace(2, 3, 5)
    x = 2*t
    y = 2*t
    y[0] = 4
    z = t
    s = t
    #t = numpy.linspace(0, 4*numpy.pi, 20)
    #cos = numpy.cos
    #sin = numpy.sin

    #x = sin(2*t)
    #y = cos(t)
    #z = cos(2*t)
    #s = 2+sin(4+t)
    
    l = points3d(2, 2, 2, 2, colormap="3", scale_factor=0.0025)
    l = points3d(2.25, 2.25, 2.25, 2.25, colormap="Spectral", scale_factor=0.001)
    l = points3d(2.5, 2.5, 2.5, 2.5, colormap="Spectral", scale_factor=0.05)
    l = points3d(2.75, 2.75, 2.75, 2.75, colormap="Spectral", scale_factor=0.006)
    l = points3d(3, 3, 3, 3, colormap="Spectral", scale_factor=0.1)
    
    
    return l

test_points3d()
show()