'''
Created on 29.09.2011

@author: axel
'''
from stats.pdistrib.sin2x_distr import sin2x
from stats.pdistrib.sin_pi_half_minus_x_distr import sin_pi_distr
import numpy as np
from math import pi as Pi
from matplotlib import pyplot as plt

x = np.linspace( 0, Pi / 2 , 10000 )
pdf = sin_pi_distr.pdf( x )
pdf[0] = 0
sinstats = sin_pi_distr.stats()
#plt.plot( x, y2 )
plt.fill( x, pdf, facecolor = 'green', alpha = 0.2 , hatch = '/' )
plt.axvline( x = sinstats[0] , linewidth = 2.0, color = 'green' , label = 'Erwartungswert' )
plt.xlabel( '$\phi$ in [$Rad$]', fontsize = 20 )
plt.ylabel( '$f_\phi (\phi)$', fontsize = 22 )
plt.title( 'Dichtefunktion PDF $sin(0,5\pi-\phi)$' , fontsize = 20 )
legend_names = ['Erwartungswert $E_\phi(\phi)$', 'Dichtefunktion $sin(0,5\pi-\phi)$']
plt.legend( legend_names, 'upper right' )
plt.xlim( 0, Pi / 2 )
plt.show()

ppf_x = np.linspace( 0, 1 , 10000 )
ppf = sin_pi_distr.ppf( ppf_x )
ppf[-1] = 0
nulls = np.zeros( len( ppf ) )
#plt.fill( x, ppf, facecolor = 'red', alpha = 0.2 , hatch = '/' )
plt.fill( ppf_x , ppf, facecolor = 'green', alpha = 0.2 , hatch = '/' )
plt.title( 'Inverse Verteilungsfunktion PPF $sin(0,5\pi-\phi)$' , fontsize = 20 )
plt.xlabel( '$F_\phi(\phi)$', fontsize = 20 )
plt.ylabel( '$\phi$ in [$Rad$]', fontsize = 22 )
plt.axhline( y = sinstats[0] , linewidth = 2.0, color = 'green' )
legend_names2 = ['Erwartungswert $E_\phi(\phi)$', 'Inverse Verteilungsfunktion $F_\phi^{-1}(\phi)$']
plt.legend( legend_names2, 'upper left' )
plt.show()
