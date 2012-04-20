from enthought.traits.api import \
    HasTraits, Float, Array, implements, Property, cached_property, Instance, \
    Int, List, Bool, HasTraits

from numpy import *

from scipy.interpolate import Rbf

from matplotlib.pyplot import *

#Title = 'Verschiebung der Elemente'
#y_label_1 = 'X'
#x = arange(0,1.01,0.01)
#y_1 = arange(0,1.01,0.01)
#
#rbf_fn = Rbf([0,.3,1],[0,.2,1], function='linear')
#
#
#color_1, color_2 = 'blue', 'red'
#
#
#fig = figure( facecolor = 'white' )
#ax1 = fig.add_subplot( 1, 1, 1 )
#
#ax1.plot( x , y_1, linewidth=1, color = 'black', label ='linear' )
#
#ax1.plot(rbf_fn(x),x,color='blue',linewidth=1,label = 'angepasst')
#ax1.plot( 0.2,0.3, 'bo')
#
#ax1.plot()
#ax1.set_xlabel( '$x_{rbf}$', fontsize=22 )
#ax1.set_ylabel( "$x'$", fontsize=22 )
#ax1.axhline(y=0.3,linewidth=1, color='b', linestyle='--')
#ax1.axvline(x=0.2,ymin=0.0, ymax=0.3,linewidth=1,color='b', linestyle='--')

Title = 'bsp_Verschiebung'
y_label_1 = 'X'
x = arange(0,1.01,0.01)
y_1 = arange(0,1.01,0.01)*4

rbf_fn = Rbf([0,.25,0.75,1],[0,0.45/2**0.5,2,4], function='linear')


color_1, color_2 = 'blue', 'red'


fig = figure( facecolor = 'white' )
ax1 = fig.add_subplot( 1, 1, 1 )

#ax1.plot( x , y_1, linewidth=1, color = 'black', label ='linear' )
ax1.plot(x,rbf_fn(x),color='blue',linewidth=1,label = 'angepasst')

x_grid = arange(0,1.01,0.25)
y_grid = rbf_fn(x_grid)

ax1.plot( [0.25,0.75],[0.45/2**0.5,2], 'bo')

#ax1.plot()
ax1.set_xlabel( '$x_{rbf}$', fontsize=22 )
ax1.set_ylabel( "$x'$", fontsize=22 )
ax1.axhline(y=0.45/2**0.5,xmax=0.25 ,linewidth=1, color='b', linestyle='--')
ax1.axhline(y=2,xmax=0.75 ,linewidth=1, color='b', linestyle='--')
ax1.axvline(x=0.25, ymax=0.45/(2**0.5)/4,linewidth=1,color='b', linestyle='--')
ax1.axvline(x=0.75, ymax=0.5,linewidth=1,color='b', linestyle='--')

ax1.set_xticks((0.25,0.75,1.0))
ax1.set_xticklabels(("$ x_{0,1} $","$ x_{0,2} $", '$1$'))
ax1.set_yticks((0.45/2**0.5,2,4.0))
ax1.set_yticklabels(("$ x'_{0,1} $","$ x'_{0,2} $",'$4$'))

for i in x_grid:
    ax1.axvline(x=i,linewidth=0.5, linestyle=':', color= 'black')

for j in y_grid:
    ax1.axhline(y=j,linewidth=0.5, linestyle=':', color= 'black')

    
for tick in ax1.xaxis.get_major_ticks():
    tick.label1.set_fontsize(15)

for tick in ax1.yaxis.get_major_ticks():
    tick.label1.set_fontsize(15)



ylim()
xlim()
#legend(loc = 5)
#for tl in ax1.get_yticklabels():
#    tl.set_color( color_1 )
#
#ax2 = ax1.twinx()
#ax2.plot( x , y_2, color = color_2 )
#ax2.set_ylabel( y_label_2, color = color_2 )
#
#for tl in ax2.get_yticklabels():
#    tl.set_color( color_2 )

fig.savefig( Title, facecolor = 'w', edgecolor = 'w',
             orientation = 'portrait' )