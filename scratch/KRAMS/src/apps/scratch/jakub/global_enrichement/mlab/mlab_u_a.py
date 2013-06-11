'''
Created on Apr 22, 2010

@author: jakub
'''

import numpy
from enthought.mayavi.mlab import surf, colorbar, show, mesh, plot3d

@show
def test_surf():
    """Test surf on regularly spaced co-ordinates like MayaVi."""
    def f(x, y):
        omega = numpy.sqrt(10.)
        sinh, cosh = numpy.sinh, numpy.cosh
        resp = numpy.zeros_like(x)
        resp[x<1.55] = 1./omega*sinh(omega*x[x<1.55])/cosh(omega*1.55)
        peak = 1./omega*sinh(omega*1.55)/cosh(omega*1.55)
        resp[x>=1.55] = 2*peak  - 1./omega*sinh(-omega*(x[x>=1.55]-3.1))/cosh(omega*1.55)
        return resp

    x, y = numpy.mgrid[0.:3.1:30j, 0.:2.1:20j]
    s = surf(x, y, f)
            #, warp_scale = 0.05)
    #cs = contour_surf(x, y, f, contour_z=0)
    return s


if __name__ == '__main__':
    test_surf()