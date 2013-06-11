'''
Created on 28.09.2011

@author: axel
'''
from stats.pdistrib.sin2x_distr import sin2x
from stats.pdistrib.sinus_distribution import sin_distr
import numpy as np
from math import pi as Pi
from matplotlib import pyplot as plt
x = np.linspace( 0, Pi / 2 , 10000 )
pdf = sin2x.pdf( x )
sin2stats = sin2x.stats()
#plt.plot( x, y2 )
plt.fill( x, pdf, facecolor = 'red', alpha = 0.2 , hatch = '\\' )
plt.axvline( x = sin2stats[0] , linewidth = 2.0, color = 'red' )
plt.xlabel( '$\phi$ in [$Rad$]', fontsize = 20 )
plt.ylabel( '$f_\phi (\phi)$', fontsize = 22 )
plt.title( 'Dichtefunktion PDF $sin(2\phi)$' , fontsize = 20 )
plt.xlim( 0, Pi / 2 )
legend_names = ['Erwartungswert $E_\phi(\phi)$', 'Dichtefunktion $sin(x)$']
plt.legend( legend_names, 'upper left' )
plt.show()

ppf_x = np.linspace( 0, 1 , 10000 )
ppf = sin2x.ppf( ppf_x )
ppf[-1] = 0
nulls = np.zeros( len( ppf ) )
#plt.fill( x, ppf, facecolor = 'red', alpha = 0.2 , hatch = '/' )
plt.fill( ppf_x , ppf, facecolor = 'red', alpha = 0.2 , hatch = '\\' )
plt.title( 'Inverse Verteilungsfunktion PPF $sin(2\phi)$' , fontsize = 20 )
plt.xlabel( '$F_\phi(\phi)$', fontsize = 20 )
plt.ylabel( '$\phi$ in [$Rad$]', fontsize = 22 )
plt.axhline( y = sin2stats[0] , linewidth = 2.0, color = 'red' )
legend_names2 = ['Erwartungswert $E_\phi(\phi)$', 'Inverse Verteilungsfunktion $F_\phi^{-1}(\phi)$']
plt.legend( legend_names2, 'upper left' )
plt.show()
