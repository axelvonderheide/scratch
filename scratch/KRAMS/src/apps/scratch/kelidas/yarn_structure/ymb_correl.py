'''
Created on Nov 19, 2010

@author: kelidas
'''
from enthought.mayavi import mlab
from enthought.mayavi.core.ui.api import MlabSceneModel, SceneEditor
from enthought.mayavi.core.ui.mayavi_scene import MayaviScene
from enthought.mayavi.tools.mlab_scene_model import MlabSceneModel
from enthought.traits.api import HasTraits, Float, Property, cached_property, \
    Event, Array, Instance, File, Int, Directory, Button, Range, Enum, \
     on_trait_change, Bool, Trait, DelegatesTo, Constant
from enthought.traits.ui.api import Item, View, Group, Handler, HSplit, VSplit, \
    HGroup, HSplit, VGroup, Tabbed, Label, Controller, ModelView
from enthought.traits.ui.api import View, Item, FileEditor, DirectoryEditor, \
    HistoryEditor
from enthought.traits.ui.menu import OKButton, CancelButton, Action, Menu, \
    MenuBar
from enthought.tvtk.pyface.scene_editor import SceneEditor
from matplotlib.figure import Figure
from numpy import loadtxt, min, array, mean, std, arange, histogram, zeros_like, \
    meshgrid, ones_like, histogram2d, c_, r_, equal, cumsum, vstack, hstack, savetxt, \
    sqrt, sum, all, zeros_like, zeros, ones, diagonal, corrcoef, linspace, sin, exp, prod
from numpy.random import random
from traits.editors.mpl_figure_editor import MPLFigureEditor
from os.path import join
from promod.simdb import SimDB
import matplotlib.pyplot as plt
import os
import re


from ymb_data import YMBData


class YMBCutCorrel( HasTraits ):

    yarn_data = Instance( YMBData, () )

    var_enum = Trait( 'contact fraction',
                      {'radius' : 'radius',
                       'shortest distance from the edge' : 'edge_distance',
                       'contact fraction' : 'cf',
                       'slack' : 'slack'} )
    # too easy
#    corr_arr = Property( Array, depends_on='var_enum, yarn_data.' )
#    @cached_property
#    def _get_corr_arr( self ):
#        corr_data = getattr( self.yarn_data, self.var_enum_ )
#        return corrcoef( corr_data, rowvar=False )

    # Continuous filaments
#    corr_arr = Property( Array, depends_on='var_enum, yarn_data.' )
#    @cached_property
#    def _get_corr_arr( self ):
#        corr_data = getattr( self.yarn_data, self.var_enum_ )
#        corr_data = corr_data[prod( corr_data >= 0, axis=1, dtype=bool )]
#        return corrcoef( corr_data, rowvar=False )

    corr_arr = Property( Array, depends_on = 'var_enum, yarn_data.' )
    @cached_property
    def _get_corr_arr( self ):
        corr_data = getattr( self.yarn_data, self.var_enum_ )
        #corr_data = corr_data[prod( corr_data >= 0, axis=1, dtype=bool )]
        import numpy.ma as ma
        corr_data = ma.masked_array( corr_data, mask = self.yarn_data.mask_arr )
        # return small differences between ma and numpy corrcoef
        #print ma.corrcoef( corr_data, rowvar=False, allow_masked=True, bias=False )
        return ma.corrcoef( corr_data, rowvar = False, allow_masked = True )

    corr_arr2 = Property()
    @cached_property
    def _get_corr_arr2( self ):
        corr = []
        for i in range( 0, self.yarn_data.n_cuts ):
            corr.append( min( 
                        corrcoef( getattr( self.yarn_data, self.var_enum_ )[:, i:-1].flatten(),
                         getattr( self.yarn_data, self.var_enum_ )[:, ( i + 1 ):].flatten() )
                        ) )
        return array( corr )

    fit_correl = Property()
    def _get_fit_correl( self ):
        x_coor = self.yarn_data.cut_x
        var_data = self.corr_arr
        x = []
        y = []
        for i in range( 0, var_data.shape[1] ):
            x.append( x_coor[i:] - x_coor[i] )
            y.append( var_data[i, ( i ):] )
        x = hstack( x )
        y = hstack( y )
        x_ls = linspace( min( x ), max( x ), 1000 )
        from scipy.optimize import leastsq
        p0 = [1., 1., 1., 1.]
        plsq = leastsq( self.residual_ls, p0, args = ( y, x ) )
        return plsq[0]

    def residual_ls( self, p, y, x ):
        err = y - self.peval( x, p )
        return err

    def peval( self, x, p ):
        return  p[0] * x ** 3 + p[1] * x ** 2 + p[2] * x + p[3]


    traits_view = View( Item( 'var_enum', label = 'Variable' ) )



class YMBCutCorrelView( HasTraits ):
    yarn_data = Instance( YMBData, () )
    correl_data = Instance( YMBCutCorrel, () )

    #n_cut = Property()
    zero = Constant( 0 )

    cut_slider = Range( 'zero', 'yarn_data.n_cuts', mode = 'slider', auto_set = False, enter_set = True, modified = True )
    vcut_slider = Range( 'zero', 'yarn_data.n_cuts', mode = 'slider', auto_set = False, enter_set = True, modified = True )

    cut_slider_on = Bool( False, modified = True )

    figure = Instance( Figure, () )

#    def _figure_default( self ):
#        figure = Figure()
#        figure.add_axes( [0.1, 0.1, 0.8, 0.8] )
#        return figure

    data_changed = Event( True )
    @on_trait_change( 'slider., +modified, correl_data.' )
    def _redraw( self ):
        # TODO: set correct ranges, fix axis range (axes.xlim)
        figure = self.figure
        figure.clear()
        var_data = self.correl_data.corr_arr
        #print mean( diagonal( self.correl_data.corr_arr, 1 ) ) 
        var_data2 = self.correl_data.corr_arr2
        #print var_data2
        id = self.cut_slider
        if self.cut_slider_on == True:
            i = self.cut_slider
            j = self.vcut_slider
            plot_data = getattr( self.yarn_data, self.correl_data.var_enum_ )
            plot_data = vstack( [plot_data[:, i], plot_data[:, j]] ).T
            # plot only non -1 values
            plot_data = plot_data[prod( plot_data >= 0, axis = 1, dtype = bool )]
            plot_data_x = plot_data[:, 0]
            plot_data_y = plot_data[:, 1]
            plot_data_corr = min( corrcoef( plot_data_x, plot_data_y ) )
            from scipy.stats import spearmanr
            plot_data_corr_spear = max( spearmanr( plot_data_x, plot_data_y ) )

            left, width = 0.1, 0.65
            bottom, height = 0.1, 0.65
            bottom_h = left_h = left + width + 0.02

            rect_scatter = [left, bottom, width, height]
            rect_histx = [left, bottom_h, width, 0.2]
            rect_histy = [left_h, bottom, 0.2, height]

            axScatter = figure.add_axes( rect_scatter )
            axHistx = figure.add_axes( rect_histx )
            axHisty = figure.add_axes( rect_histy )
            axScatter.clear()
            axHistx.clear()
            axHisty.clear()

            from matplotlib.ticker import NullFormatter
            axHistx.xaxis.set_major_formatter( NullFormatter() )
            axHisty.yaxis.set_major_formatter( NullFormatter() )

            axScatter.scatter( plot_data_x,
                               plot_data_y )

            #binwidth = 0.25
            #xymax = max( [max( abs( self.yarn_data.cf[:, j] ) ), max( abs( self.yarn_data.cf[:, i] ) )] )
            #lim = ( int( xymax / binwidth ) + 1 ) * binwidth

            #axScatter.set_xlim( ( -lim, lim ) )
            #axScatter.set_ylim( ( -lim, lim ) )

            #bins = arange( -lim, lim + binwidth, binwidth )
            axHistx.hist( plot_data_x, bins = 40 )
            axHisty.hist( plot_data_y, bins = 40, orientation = 'horizontal' )
            axHistx.set_xlim( axScatter.get_xlim() )
            axHisty.set_ylim( axScatter.get_ylim() )

            axScatter.set_xlabel( '$\mathrm{cut\, %i}$' % self.cut_slider, fontsize = 16 )
            axScatter.set_ylabel( '$\mathrm{cut\, %i}$' % self.vcut_slider, fontsize = 16 )
            axScatter.text( axScatter.get_xlim()[0], axScatter.get_ylim()[0],
                             'actual set correlation %.3f (Pearson), %.3f (Spearman)' % ( plot_data_corr, plot_data_corr_spear ), color = 'r' )

        if self.cut_slider_on == False:
            figure.add_axes( [0.1, 0.1, 0.8, 0.8] )
            axes = figure.axes[0]
            axes.clear()
            x_coor = self.yarn_data.cut_x
            axes.grid()
            for i in range( 0, var_data.shape[1] ):
                axes.plot( x_coor[i:] - x_coor[i] ,
                           var_data[i, ( i ):], '-x' )
            axes.plot( x_coor, self.correl_data.peval( x_coor, self.correl_data.fit_correl ), 'b', linewidth = 3 )
            axes.set_xlabel( '$\mathrm{x}\, [\mu\mathrm{m}]$', fontsize = 16 )
            axes.set_ylabel( '$\mathrm{correlation}$', fontsize = 16 )
            axes.set_ylim( -1, 1 )

#        if self.cut_slider_on == False:
#            for i in range( 0, var_data.shape[1] ):
#                axes.plot( arange( i, len( var_data[i] ) ),
#                           var_data[i, ( i ):] )
        self.data_changed = True


    traits_view = View( Group( Item( 'correl_data', show_label = False, style = 'custom' ),
                              HGroup( 
                       Item( 'cut_slider_on', label = 'Scatter' ),
                       Item( 'cut_slider', show_label = False, springy = True, enabled_when = 'cut_slider_on == True' ),
                       Item( 'vcut_slider', show_label = False, springy = True, enabled_when = 'cut_slider_on == True' ),
                       ),
                       Group( Item( 'figure', style = 'custom',
                              editor = MPLFigureEditor(),
                              show_label = False )
                              , id = 'figure.view' ),
                       ),
                       resizable = True,
                        )



if __name__ == '__main__':
    corr = YMBCutCorrelView( correl_data = YMBCutCorrel() )
    corr.configure_traits()

    from matplotlib.colors import LinearSegmentedColormap, Colormap

    from numpy import mgrid, max
    mg = mgrid[0:100:36j, 0:100:36j]
    d = corr.yarn_data.cf
    d2 = d[:, 0:2][prod( d[:, 0:2] >= 0, axis = 1, dtype = bool )]
    a, b, c = histogram2d( d2[:, 0], d2[:, 1], bins = 36 )

    # my own colormap
    cdict = {
    'red':   ( ( 0.0, 1.0, 1.0 ),
                ( 1 / max( a ) , 1.0, 1.0 ),
                 ( 1 / max( a ), 0.0, 0.0 ),
                  ( 1.0, 0.0, 0.0 ) ),
    'green': ( ( 0.0, 0.0, 0.0 ),
                ( 1 / max( a ) , 0.0, 0.0 ),
                 ( 1 / max( a ), 1.0, 1.0 ),
                  ( 1.0, 0.0, 0.0 ) ),
    'blue':  ( ( 0.0, 0.0, 0.0 ),
                ( 1 / max( a ) , 0.0, 0.0 ),
                 ( 1 / max( a ), 0.0, 0.0 ),
                  ( 1.0, 1.0, 0.5 ) ) }

    my_cmap_lin = LinearSegmentedColormap( 'my_colormap_lin', cdict, 256 )
    plt.title( 'Contact fraction scatter -- cut0 x cut1' )
    scat = plt.scatter( mg[0], mg[1], c = a, s = 40, cmap = my_cmap_lin )
    plt.colorbar( scat )
    plt.show()

