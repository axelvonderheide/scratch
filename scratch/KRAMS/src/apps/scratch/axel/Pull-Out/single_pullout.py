'''
Created on 17.09.2011

@author: axel
'''
from quaducom.resp_func.cb_short_fiber import CBShortFiber
from numpy import linspace, argmax, zeros
from matplotlib import pyplot as plt
from math import pi

Cbsf = CBShortFiber()
w = linspace( 0, 1, 10000 )
all_resp = Cbsf( w, 2.3, 4.5, 2 * 0.075, 200000, 4.5, pi / 2, 1, 0, 0, 1e15 )
plt.ylabel( 'Kraft $P$ in [$N$]' , fontsize = 22 )
plt.xlabel( 'Verschiebung $w$ in [$mm$]' , fontsize = 22 )
plt.ylim( 0, 38 )
plt.xlim( 0, 10 )
plt.title( 'Phasen des Faserauszugs' , fontsize = 23 )


divide_index = argmax( all_resp )

print divide_index
all_resp[divide_index ] = 0
null = zeros( len( w[divide_index + 1:-1] ) )
plt.fill( w[0:divide_index + 1] , all_resp[0: divide_index + 1], facecolor = 'blue', alpha = 0.35 )
plt.fill_between( w[divide_index + 1:-1] , all_resp[divide_index + 1:-1], null, facecolor = 'red', alpha = 0.35 )
plt.fill( [-1, -2] , [1, 2], facecolor = 'red', alpha = 0.35 )
legend_names2 = [  'Faserauszug', 'Faseraktivierung']
plt.legend( legend_names2, 'upper left' )
plt.axvline( x = w[divide_index], linewidth = 5, color = 'black' , ls = '--' )
plt.plot( w, all_resp, 'k' )
plt.xlim( 0, 1 )
plt.ylim( 0, 30 )
plt.show()
