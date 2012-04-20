from numpy import *
from matplotlib import pyplot

#def bar_plot( x, y, filename = 'bla', title = 'Title', xlabel = 'xlabel', ylabel = 'ylavel', width = 0.1 ):
#            print x
#            print y
#            fig = pyplot.figure( facecolor = "white" , figsize = [10, 5] )
#            ax1 = fig.add_subplot( 1, 1, 1 )
#            ax1.bar( x , y , width = width, align = 'center', color = 'blue' )
#            ax1.set_xlabel( xlabel, fontsize = 22 )
#            ax1.set_ylabel( ylabel, fontsize = 22 )
#            ax1.set_title( 'TITLE' )
#            fig.savefig( filename, orientation = 'portrait', bbox = '' )
#            pyplot.show()
#            pyplot.clf()
#
#def multi_bar_plot( x, y_1, y_2, y_3, filename = 'bla', title = 'Title', xlabel = 'xlabel',
#                   y_1_label = 'ylabel1', y_2_label = 'ylabel2', y_3_label = 'ylabel3', width = 0.1 ):
#
#            fig = pyplot.figure( facecolor = "white" , figsize = [10, 5] )
#            ax1 = fig.add_subplot( 3, 1, 1 )
#            ax1.bar( x , y_1 , width = width, align = 'center', color = 'blue' )
#            ax1.set_ylabel( y_1_label, fontsize = 22 )
#            ax1.set_title( 'TITLE' )
#            ax2.subplot( 312, sharex = ax1 )
#            ax2.bar( x , y_2, width = width, align = 'center', color = 'blue' )
#            ax2.set_ylabel( y_2_label, fontsize = 22 )
#            ax3.subplot( 313, sharex = ax1 )
#            ax3.bar( x , y_3 , width = width, align = 'center', color = 'blue' )
#            ax3.set_ylabel( y_3_label, fontsize = 22 )
#            show()
##
#def interaction_plot():
#
#        fig = pyplot.figure( facecolor = "white", figsize = [10, 10] )
#        ax1 = fig.add_subplot( 1, 1, 1 )
#        x = arange( 0, 1.01, 0.01 )
#        y = ( 1 - x ** 1.5 ) ** ( 1 / 1.5 )
#
#        ax1.set_xlabel( '$V_{Ed}/V_{Rd}$' , fontsize = 22 )
#        ax1.set_ylabel( '$N_{Ed}/N_{Rd}$', fontsize = 22 )
#        ax1.plot( x , y, '--', color = 'black', label = 'quadratische Interaktion' )
##        ax1.plot( x , 1 - x, '--', color = 'black', label = 'lineare Interaktion' )
#
#        ax1.set_xlim( 0, 1.1 )
#        ax1.set_ylim( 0, 1.1 )
#        ax1.legend()
#        pyplot.show()
#        pyplot.clf()
#
#a = arange(0,21,1).reshape(7,3)
#print a
#b=array([0,9,15])

#idx = []
#for bs in b:
#    if (abs(a[:,0] - bs)  <= 0.0001).any()==True:
#    print list(where(abs(a[:,0] - bs) <= 0.0001)[0])
#    idx = idx + list(where(a[:,0] == bs)[0])
#    print a[idx,:]

#print idx 


#a = arange( 0, 10, 1 )
#b = a * 2
#print a
#print b
#multi_bar_plot( a, b, b, b, xlabel = '$X$ [m]' )
#interaction_plot()

#x = arange(0,48,1).reshape(6,2,4)
#print x
#print x[:,:,0]
#print argmax(x[:,:,0],axis=0).shape
#print x[argmax(x[:,:,0],axis=1)].shape

a = 'bla'
# % ( xkey, )
b = '$ %s $' % ( a )
print b
