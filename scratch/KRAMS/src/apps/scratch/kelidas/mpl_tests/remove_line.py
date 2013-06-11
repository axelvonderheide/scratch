'''
Created on Jan 24, 2011

@author: kelidas
'''

'''
    removing one data series from an already existing plot
'''

import matplotlib.pyplot as plt
from numpy import arange, array

x = arange( 10 )

fig0 = plt.figure( 0 )

ax0 = fig0.add_subplot( 111 )
ax0.plot( x )

ax0.plot( x + 10 )

ax0.plot( x + 20 )



print ax0.lines


fig1 = plt.figure( 1 )

ax1 = fig1.add_subplot( 111 )
#del ax0.lines[1]
ax1.lines = [ax0.lines[0], ax0.lines[2]]

plt.show()
