'''
Created on 28.09.2011

@author: axel
'''
from stats.pdistrib.sin2x_distr import sin2x
from stats.pdistrib.sinus_distribution import sin_distr
from scipy.stats import norm
import numpy as np
from math import pi as Pi
from matplotlib import pyplot as plt
x = np.linspace( norm( 763.94, 27.54 ).ppf( 0.00001 ), norm( 763.94, 27.54 ).ppf( 0.99999 ) , 10000 )
pdf = norm( 763.94, 27.54 ).pdf( x )
#sin2stats = sin2x.stats()
#plt.plot( x, y2 )
plt.fill( x, pdf, facecolor = 'green', alpha = 0.4 , hatch = '|' )
#plt.axvline( x = sin2stats[0] , linewidth = 2.0, color = 'red' )
plt.xlabel( 'Anzahl der Fasern $n$ ', fontsize = 24 )
plt.ylabel( '$f_n (n)$', fontsize = 24 )
plt.title( 'Dichtefunktion PDF Normalverteilung' , fontsize = 20 )
plt.xlim( 663.94, 863.94 )
legend_names = ['$E_x(x)=763.94, \sigma_x(x)=27.54$']
plt.legend( legend_names, 'upper left' )
plt.show()

ppf_x = np.linspace( 0, 1 , 1000 )
ppf = norm( 763.94, 27.54 ).ppf( ppf_x )
ppf[0] = 0
ppf[-1] = 0
nulls = np.zeros( len( ppf ) )
#plt.fill( x, ppf, facecolor = 'red', alpha = 0.2 , hatch = '/' )
plt.fill( ppf_x , ppf, facecolor = 'green', alpha = 0.4 , hatch = '|' )
plt.title( 'Inverse Verteilungsfunktion PPF Normalverteilung' , fontsize = 24 )
plt.xlabel( '$F_n(n)$', fontsize = 20 )
plt.ylabel( 'Anzahl der Fasern $n$ ', fontsize = 24 )
plt.axhline( y = 763.94 , linewidth = 2.0, color = 'green' )
legend_names2 = ['Erwartungswert $E_n(n)$', 'Inverse Verteilungsfunktion $F_n^{-1}(n)$']
plt.legend( legend_names2, 'upper left' )
plt.ylim( 663.94, 863.94 )
plt.show()

