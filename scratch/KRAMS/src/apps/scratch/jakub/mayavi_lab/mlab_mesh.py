'''
Created on Sep 2, 2009

@author: jakub
'''
import numpy
from enthought.mayavi.mlab import surf, colorbar, show, mesh, plot3d

@show
def mmesh():
    """A very pretty picture of spherical harmonics translated from
    the octaviz example."""
    pi = numpy.pi
    cos = numpy.cos
    sin = numpy.sin
    dphi, dtheta = pi/250.0, pi/250.0
    [phi,theta] = numpy.mgrid[0:pi+dphi*1.5:dphi,0:2*pi+dtheta*1.5:dtheta]
    m0 = 4; m1 = 3; m2 = 2; m3 = 3; m4 = 6; m5 = 2; m6 = 6; m7 = 4;
    r = sin(m0*phi)**m1 + cos(m2*phi)**m3 + sin(m4*theta)**m5 + cos(m6*theta)**m7
    x = r*sin(phi)*cos(theta)
    y = r*cos(phi)
    z = r*sin(phi)*sin(theta);

    return mesh(x, y, z, colormap="bone")

@show
def gsurf(): 
    """Test contour_surf on regularly spaced co-ordinates like MayaVi."""
    def f(x, y):
        return x**2+y**2 - 25.

    x, y = numpy.mgrid[-10.:10.:0.25, -10.:10.0:0.25]
    s = surf(x, y, f, warp_scale = 0.)
    colorbar(s,orientation='horizontal')
    return s

@show
def line():
    """Generates a pretty set of lines."""
    x =y = z = numpy.arange(5)

    l = plot3d(x, y, z, 
               #tube_radius=0.025, 
               color=(0.862745,0.0784314,0.235294))




if __name__ == '__main__':
    gsurf()
