from numpy import abs, sin, cos, mgrid, pi
from enthought.mayavi.scripts import mayavi2
from enthought.mayavi import mlab

if __name__ == '__main__':
        
    mlab.figure(fgcolor=(0, 0, 0), bgcolor=(0.8, 0.1, 0.1))
    u, v = mgrid[0.:2*pi:0.2, 0.:2*pi:0.2]
    r = 0.3
    
    X = r* cos(u)* sin(v)
    Y = r* sin(u)* sin(v)
    Z = r* cos(v)
    S = abs(X)*abs(Y)
    print u,v
    
    mlab.mesh(X, Y, Z, scalars=S, colormap='gist_earth', )
    
    XX = 0.7*r* cos(u)* sin(v)
    YY = 0.7 + 0.7*r* sin(u)* sin(v)
    ZZ = 0.7*r* cos(v)
    SS = abs(XX)*cos(YY)
    
    mlab.mesh(XX, YY, ZZ, scalars=SS, colormap='Spectral', )
    
    mlab.view(.0, -5.0, 2.5)
    mlab.show()
