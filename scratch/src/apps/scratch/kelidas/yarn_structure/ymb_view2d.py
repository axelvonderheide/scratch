'''
Created on Nov 20, 2010

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
    HistoryEditor, ShellEditor
from enthought.traits.ui.menu import OKButton, CancelButton, Action, Menu, \
    MenuBar
from enthought.tvtk.pyface.scene_editor import SceneEditor
from matplotlib.figure import Figure
from numpy import loadtxt, min, array, mean, std, arange, histogram, zeros_like, \
    meshgrid, ones_like, histogram2d, c_, r_, equal, cumsum, vstack, hstack, savetxt, \
    sqrt, sum, all, zeros_like, zeros, ones, max, ceil
from numpy.random import random
from traits.editors.mpl_figure_editor import MPLFigureEditor
from os.path import join
from promod.simdb import SimDB
import matplotlib.pyplot as plt
from enthought.pyface.api import ApplicationWindow, GUI, PythonShell
import os
import re
from scipy.stats import corrcoef
from matplotlib.colors import LinearSegmentedColormap, Colormap

from ymb_data import YMBData


# my own colormap, similar to ATENA -- deafault = not used
cdict = {
    'red':   ( ( 0.0, 1.0, 1.0 ), ( 0.5, 0.0, 0.0 ), ( 1.0, 0.0, 0.0 ) ),
    'green': ( ( 0.0, 0.0, 0.0 ), ( 0.5, 1.0, 1.0 ), ( 1.0, 0.0, 0.0 ) ),
    'blue':  ( ( 0.0, 0.0, 0.0 ), ( 0.5, 0.0, 0.0 ), ( 1.0, 1.0, 0.5 ) )
    }
my_cmap_lin = LinearSegmentedColormap( 'my_colormap_lin', cdict, 256 )


class YMBView2D( HasTraits ):
    data = Instance( YMBData, () )

    zero = Constant( 0 )
    n_cuts = DelegatesTo( 'data' )
    n_filaments = DelegatesTo( 'data' )

    var_enum = Trait( 'radius',
                      {'radius' : 'radius',
                       'area':'area',
                       'shortest distance from the edge' : 'edge_distance',
                       'contact fraction' : 'cf'
                        }, modified=True )

    cut_slider = Range( 'zero', 'n_cuts', mode='slider', auto_set=False,
                         enter_set=True, modified=True )

    circle_diameter = Float( 20, enter_set=True, auto_set=False, modified=True )

    underlay = Bool( False, modified=True )

    variable = Property( Array, depends_on='var_enum, data.' )
    @cached_property
    def _get_variable( self ):
        return getattr( self.data, self.var_enum_ )

    figure = Instance( Figure, () )

    def _figure_default( self ):
        figure = Figure()
        figure.add_axes( [0.1, 0.1, 0.8, 0.8] )
        return figure

    data_changed = Event( True )
    @on_trait_change( ' +modified' )
    def _redraw( self ):
        # TODO: set correct ranges, fix axis range (axes.xlim)
        self.figure.clear()
        self.figure.add_axes( [0.1, 0.1, 0.8, 0.8] )
        figure = self.figure
        axes = figure.axes[0]
        axes.clear()
        x_arr, y_arr, z_arr, r_arr, d_arr, cf_arr = self.data.cut_data
        y_raw_arr = self.data.cut_raw_data[0]
        z_raw_arr = self.data.cut_raw_data[1]
        offset = hstack( [0, self.data.cut_raw_data[5]] )
        print self.data.cut_raw_data[5]
        print offset
        scalar_arr = self.variable
        mask = y_arr[:, self.cut_slider] > -1

        scat_raw = axes.scatter( y_raw_arr[offset[self.cut_slider]:offset[self.cut_slider + 1]],
                                  z_raw_arr[offset[self.cut_slider]:offset[self.cut_slider + 1]],
                                  s=self.circle_diameter, color='k', marker='x', label='identified filament in cut' )
        #savetxt( 'cut11_raw.txt', vstack( [y_raw_list[self.cut_slider],
        #                          z_raw_list[self.cut_slider], ] ).T )
        #savetxt( 'cut11.txt', vstack( [y_arr[:, self.cut_slider][mask], z_arr[:, self.cut_slider][mask], scalar_arr[:, self.cut_slider][mask]] ).T )
        scat = axes.scatter( y_arr[:, self.cut_slider][mask], z_arr[:, self.cut_slider][mask],
                      s=self.circle_diameter, c=scalar_arr[:, self.cut_slider][mask], cmap=my_cmap_lin, label='connected filaments' )
        axes.set_xlabel( '$y\, [\mu\mathrm{m}]$', fontsize=16 )
        axes.set_ylabel( '$z\, [\mu\mathrm{m}]$', fontsize=16 )
        #axes.set_xlim( [0, ceil( max( y_arr ) / 100 ) * 100] )
        #axes.set_ylim( [0, ceil( max( z_arr ) / 100 ) * 100] )
        axes.legend()
        figure.colorbar( scat )

        if self.underlay == True:
            axes.text( axes.get_xlim()[0], axes.get_ylim()[0],
                        'That\'s all at this moment :-)', color='red', fontsize=20 )
        self.data_changed = True


    traits_view = View( Group( Item( 'var_enum' ),
                       Item( 'cut_slider', springy=True ),
                       Item( 'circle_diameter', springy=True ),
                       Item( 'underlay', springy=True ),
                        ),
                        Item( 'figure', style='custom',
                              editor=MPLFigureEditor(),
                              show_label=False ),
                       resizable=True,
                        )




if __name__ == '__main__':
    cut = YMBView2D()
    cut.configure_traits()
