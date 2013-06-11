
import numpy
from enthought.mayavi.mlab import *
from enthought.mayavi import mlab
import math

@mlab.show
def test_surf():
    """Test surf on regularly spaced co-ordinates like MayaVi."""

    R = 0.0006
    Alpha_min = 0.0
    Alpha_max = math.pi / 2.0
    Alpha = 2.0 * math.pi
    n = 3j
    
    def f(r, alpha):
        return math.sqrt( r ) # * math.cos( alpha )

    r, alpha = numpy.mgrid[0:R:5j, Alpha_min : Alpha_max : n]

    x = r * numpy.cos( alpha )
    y = r * numpy.sin( alpha )

    sfunc = numpy.frompyfunc( f, 2, 1 )
    z = numpy.array( sfunc( r, alpha ), dtype = 'float_' )

    print 'x',x
    print 'y',y
    print 'z',z 


    s = mesh(x, y, z, colormap="bone" )
    #cs = contour_surf(x, y, f, contour_z=0)
    return s

s = test_surf()