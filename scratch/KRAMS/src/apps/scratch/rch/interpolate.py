from numpy import *
from pylab import *
#from scikits.delaunay import *
from numpy import ma

def plot_nodes(tri):
    for nodes in tri.triangle_nodes:
        fill(x[nodes],y[nodes],'b')
    show()

def plot_data(xi,yi,zi):
    zim = ma.masked_where(isnan(zi),zi)
    figure(figsize=(8,8))
    pcolor(xi,yi,zim,shading='interp',cmap=cm.gray)
    contour(xi,yi,zim,cmap=cm.jet)
    show()

if __name__ == '__main__':
    N = 10000
    aspect = 1.0

    # Data
    x = randn(N)/aspect
    y = randn(N)
    z = rand(N)

    # Grid
    xi, yi = mgrid[-5:5:100j,-5:5:100j]

    # triangulate data
    tri = Triangulation(x,y)

    # interpolate data
    interp = tri.nn_interpolator(z)

    zi = interp(xi,yi)
    # or, all in one line
    #    zi = Triangulation(x,y).nn_interpolator(z)(xi,yi)

    plot_data(xi,yi,zi)