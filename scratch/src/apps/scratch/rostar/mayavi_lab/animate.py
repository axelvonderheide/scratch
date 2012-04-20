from numpy import array, mgrid, vstack, arange, transpose, zeros
from enthought.mayavi.mlab import show, surf, view, roll, axes, xlabel


x, y = mgrid[0:3:1,0:3:1]
s = surf(x, y, zeros(9).reshape(3,3), warp_scale = 0.0, vmin = 0, vmax = 3)
view(0,0,20)
for i,j in transpose(vstack((arange(100),arange(100)))):
    s.mlab_source.x = x*(1+0.01*(i+1))
    s.mlab_source.y = y*(1+0.003*(j+1))
    s.mlab_source.scalars = (x*(1+0.01*i)-x)+(y*(1+0.01*j)-y)
axes()
show()


