'''
Created on 28.09.2011

@author: axel
'''

from math import e
from matplotlib import pyplot as plt
import numpy as np
lcorr = 6
def acor( dx, lcorr ):
    '''autocorrelation function'''
    return e ** ( -( dx / lcorr ) ** 2 )

x = np.linspace( 0, 16, 1000 )
y = acor( x, lcorr )
y[0] = 0

plt.fill( x, y, color = 'blue' , alpha = .2 )
plt.xlabel( 'Abstand zu $x$ in [mm]' , fontsize = 18 )
plt.ylabel( '$F_lcorr(x)$ Korrelationsfaktor [-]' , fontsize = 18 )
name = ['Autokorellationsfunktion']

plt.legend( name )
plt.title( 'Autokorrelation mit lcorr=6', fontsize = 20 )
plt.plot( x, y , 'k' )
plt.show()
