
from mathkit.mfn.mfn_line.mfn_line import MFnLineArray
from numpy import cos, sin, linspace

x = linspace( 0., 1., 100 )
y = 1.0 + 0.5 * sin( 9.0 * x ) * cos( 10 * x )

print x
print y

mfn =  MFnLineArray( xdata = x, ydata = y )

mfn.configure_traits()