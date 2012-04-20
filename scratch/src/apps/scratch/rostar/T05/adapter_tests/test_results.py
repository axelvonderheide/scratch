'''
Created on Aug 8, 2011

@author: rostar
'''

from enthought.traits.api import \
    HasTraits, Str, List, Tuple

import numpy as np

from scipy.interpolate import RectBivariateSpline as spline

from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt

class AdapterTests( HasTraits ):

    title = Str
    label1 = Str
    label2 = Str

    name1 = Str
    name2 = Str
    name3 = Str

    xlim = Tuple
    ylim = Tuple
    zlim = Tuple

    values1 = List
    values2 = List
    parameter = List

    results = List

    def param_plot( self ):
        fig = plt.figure()
        ax = fig.gca( projection = '3d' )
        values = self.results
        p1 = self.values1
        p2 = self.values2
        n1 = self.name1
        n2 = self.name2
        n3 = self.name3

        x = np.linspace( p1[0], p1[-1], 50 )
        y = np.linspace( p2[0], p2[-1], 50 )

        X = np.meshgrid( x, x )[0]
        Y = np.meshgrid( y, y )[1]

        val = np.array( values[0:9] ).reshape( 3, 3 )
        spl = spline( p1, p2, val, kx = 2, ky = 2 )
        Z = spl( x, y ) / 0.89

        val2 = np.array( values[9:] ).reshape( 3, 3 )
        spl2 = spline( p1, p2, val2, kx = 2, ky = 2 )
        ZZ = spl2( x, y ) / 0.89

        ax.plot_wireframe( X, Y, Z, rstride = 5, cstride = 5, alpha = 1.,
                         color = 'black', lw = 2 )

        i = 0
        for V in p2:
            for vorsp in p1:
                print vorsp, V, values[i]
                ax.plot( [vorsp, vorsp], [V, V], values[i] / 0.89, 'ko', lw = 5 )
                i += 1
        # glas 1200tex 5cN: ax.plot( [100, 100], [30, 30], np.array( [588.32, 588.32] ) / 0.89, 'ro', label = '%.1f MPa' % ( 588.32 / 0.89 ) )
        ax.plot( [100, 100], [30, 30], np.array( [588.32, 588.32] ) / 0.89, 'ro', label = '%.1f MPa' % ( 588.32 / 0.89 ) )

        ax.plot( [0, 0], [-1e15, -1e15], color = 'black', label = self.label1, lw = 2 )


#        ax.plot_wireframe( X, Y, ZZ, rstride = 5, cstride = 5, alpha = 1.,
#                         color = 'black', lw = 2 )
#
#        i = 0
#        for V in p2:
#            for vorsp in p1:
#                print vorsp, V, values[9 + i]
#                ax.plot( [vorsp, vorsp], [V, V], values[9 + i] / 0.89, 'ko' )
#                i += 1
#        ax.plot( [50, 50], [30, 30], np.array( [594.15, 594.15] ) / 0.89, 'ro', label = '%.1f MPa' % ( 594.15 / 0.89 ) )
#
#        ax.plot( [0, 0], [-1e15, -1e15], color = 'black', label = self.label2, lw = 2 )


        ax.legend( loc = 'upper left' )
        ax.set_title( self.title )
        ax.set_xlabel( n1 )
        ax.set_xlim3d( self.xlim[0], self.xlim[1] )
        ax.set_ylabel( n2 )
        ax.set_ylim3d( self.ylim[0], self.ylim[1] )
        ax.set_zlabel( n3 )
        ax.set_zlim3d( self.zlim[0], self.zlim[1] )

        plt.show()

    lengths = List
    force = List
    cov = List

    def se_plot( self ):
        l = self.lengths
        strength = np.array( self.force ) / 0.89
        stdev = np.array( self.cov ) / 100. * strength
        plt.plot( l, strength, color = 'black', lw = 2 )
        plt.errorbar( elinewidth = 2, x = l, y = strength,
                   yerr = stdev, color = 'black', fmt = 'o' )
        plt.ylim( ymin = 0.8 * np.min( strength ), ymax = 1.1 * np.max( strength ) )
        plt.xlim( xmax = 600 )
        plt.xlabel( 'Laenge [mm]' )
        plt.ylabel( 'Festigkeit [MPa]' )
        plt.grid()
        plt.title( self.title )
        plt.show()

if __name__ == '__main__':
    at = AdapterTests()

    def glas2400params():
        at.title = 'AR-Glas 2400 tex'
        at.label1 = '$P_0 = 100\%$'
        at.label2 = '$P_0 = 60\%$'

        at.name1 = 'Vorspannung'
        at.name2 = 'V-Kraft'
        at.name3 = 'Festigkeit [MPa]'

        at.xlim = ( 0, 15 )
        at.ylim = ( 20, 60 )
        at.zlim = ( 950, 1150 )

        at.parameter = [100., 60.]
        at.values1 = [1., 7., 15.]
        at.values2 = [25., 40., 55.]
        at.results = [920.29, 965.93, 961.84, 957.62, 978.63,
                  998.54, 940.93, 963.20, 963.45, 900.66,
                  864.42, 911.20, 918.08, 970.91, 971.69,
                  905.19, 956.57, 970.19]


        at.param_plot()

    def glas2400se():
        at.title = 'AR-Glas 2400 tex'
        at.lengths = [50., 70, 100, 150, 200, 250, 350, 500]
        at.force = [920.04, 912.50, 913.92, 911.31, 922.06, 898.67, 893.16, 878.78]
        at.cov = [4.30, 4.44, 3.14, 2.89, 4.76, 4.39, 5.17, 4.40]
        at.se_plot()

    def glas1200params():

        at.title = 'AR-Glas 1200 tex'
        at.label1 = 'Vorspannung = 5cN/tex'
        at.label2 = 'Vorspannung = 15cN/tex'

        at.name1 = 'Anpresskraft'
        at.name2 = 'V-Kraft'
        at.name3 = 'Festigkeit [MPa]'

        at.xlim = ( 50, 100 )
        at.ylim = ( 15, 50 )
        at.zlim = ( 600, 680 )

        at.parameter = [5., 15.]
        at.values1 = [50., 80., 100.]
        at.values2 = [15., 30., 50.]
        at.results = [562.92, 546.85, 571.11, 573.47,
                      574.42, 588.32, 572.03, 574.25,
                      571.14, 547.39, 562.03, 550.87,
                      594.15, 582.88, 585.40, 579.56,
                      577.29, 583.20]

        at.param_plot()

    def glas1200se():
        at.title = 'AR-Glas 1200 tex'
        at.lengths = [50., 70, 100, 150, 200, 250, 350, 500]
        at.force = [579.36, 570.41, 578.92, 572.17, 580.79, 567.55, 574.63, 550.96]
        at.cov = [2.74, 2.78, 4.31, 3.57, 4.47, 4.20, 4.63, 5.32]
        at.se_plot()

    def carbon1600params():

        at.title = 'Carbon 1600 tex'
        at.label1 = ''
        at.label2 = ''

        at.name1 = ''
        at.name2 = ''
        at.name3 = 'Festigkeit [MPa]'

        at.xlim = ()
        at.ylim = ()
        at.zlim = ()

        at.parameter = []
        at.values1 = []
        at.values2 = []
        at.results = []

        at.param_plot()

    def carbon1600se():
        at.title = 'Carbon 1600 tex'
        at.lengths = [50., 70, 100, 150, 200, 250, 350, 500]
        at.force = [1720.37, 1737.90, 1714.25, 1526.31, 1477.59, 1381.55, 1269.31, 1246.36]
        at.cov = [5.16, 3.13, 4.37, 7.14, 8.69, 9.23, 7.12, 10.81]
        at.se_plot()

    def carbon800params():

        at.title = 'Carbon 1600 tex'
        at.label1 = ''
        at.label2 = ''

        at.name1 = ''
        at.name2 = ''
        at.name3 = 'Festigkeit [MPa]'

        at.xlim = ()
        at.ylim = ()
        at.zlim = ()

        at.parameter = []
        at.values1 = []
        at.values2 = []
        at.results = []

        at.param_plot()

    def carbon800se():
        at.title = 'Carbon 800 tex'
        at.lengths = [50., 70, 100, 150, 200, 250, 350, 500]
        at.force = [990.08, 936.61, 880.82, 846.06, 819.32, 719.08, 708.00, 656.11]
        at.cov = [3.44, 5.07, 5.75, 9.94, 8.98, 8.30, 7.28, 9.16]
        at.se_plot()

    #glas2400params()
    #glas2400se()
    #glas1200params()
    #glas1200se()

    #carbon1600params()
    #carbon1600se()
    #carbon800params()
    carbon800se()
