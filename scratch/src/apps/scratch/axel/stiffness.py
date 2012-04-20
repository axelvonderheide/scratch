from numpy import array, sqrt
from numpy import array, mean, var, sqrt, linspace
from matplotlib import pyplot as p
from scipy.optimize import curve_fit

f = open( '/tensile_tests/Masche.raw', 'r' )
i = 0
list_x = []
list_y = []

for line in f:
    i += 1
    if i > 41:
        #X-Data
        tmp1 = line[line.index( ';' ) + 1:]
        tmp1 = tmp1.replace( ',', '.' )
        strain = tmp1[:tmp1.index( ';' )]
        strain = strain.replace( ',', '.' )
        list_x.append( float( strain ) )
        #Y-Data
        tmp2 = tmp1[tmp1.index( ';' ) + 1:]
        tmp2 = tmp2[tmp2.index( ';' ) + 1:]
        strength = float( tmp2[:tmp2.index( ';' )] )
        list_y.append( strength )

if list_x[0] < 0:
    list_x += sqrt( list_x[0] ** 2 )
else: list_x -= sqrt( list_x[0] ** 2 )
shortx = list_x[170:300]
shorty = list_y[170:300]
def fitstiffness( x, a, b ):
    return a * x + b
a, varf = curve_fit( fitstiffness, shortx, shorty )
print a
#p.plot( list_x, list_y, 'b' )
x_array = linspace( 0, 5, 50 )
y = a[0] * x_array + a[1]
p.plot( shortx, shorty , 'b' )
p.plot( x_array, y, 'y' )
p.show()
f.close()
