
from mathkit.mfn.mfn_line.mfn_line import MFnLineArray
from numpy import linspace, frompyfunc, array
from math import sin, cos

sigma_max = 0.4

xdata = linspace( 0., 1., 10000 )
fnydata = lambda x:  1.0+0.0001*sin(19.0*x)*cos(200*x) * sigma_max * 10000

fnydata_vec = frompyfunc( fnydata, 1, 1 )
ydata = array( fnydata_vec( xdata ), dtype = float )

#print ydata

mfn = MFnLineArray( xdata = xdata, ydata = ydata )
mfn.configure_traits()
