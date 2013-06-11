'''
Created on Dec 13, 2010

@author: kelidas
'''
from enthought.traits.api import HasTraits, Float, Property, cached_property, \
    Event, Array, Instance, File, Int, Directory, Button, Range, Enum, \
     on_trait_change, Bool, Trait, HasPrivateTraits, Constant, List, Tuple, \
     DelegatesTo, Interface, implements
from enthought.traits.ui.api import Item, View, Group, Handler, HSplit, VSplit, \
    HGroup, HSplit, VGroup, Tabbed, Label, Controller, ModelView, \
    FileEditor, DirectoryEditor, \
    HistoryEditor, RangeEditor
from enthought.traits.ui.menu import OKButton, CancelButton, Action, Menu, \
    MenuBar
from enthought.tvtk.pyface.scene_editor import SceneEditor
from matplotlib.figure import Figure
from numpy import loadtxt, min, array, mean, std, arange, histogram, zeros_like, \
    meshgrid, ones_like, histogram2d, c_, r_, equal, cumsum, vstack, hstack, savetxt, \
    sqrt, sum, all, zeros_like, zeros, ones, diff, where, unique, isnan, pi, invert, \
    histogram, trapz, linspace
from numpy.random import random
from os.path import join
from promod.simdb import SimDB
from traits.editors.mpl_figure_editor import MPLFigureEditor
import matplotlib.pyplot as plt
import numpy.ma as ma
import os
import re
from enthought.traits.trait_types import DelegatesTo
from mathkit.mfn.mfn_line.mfn_line import MFnLineArray

from ymb_data import YMBSlider
from ymb_data import YMBData
from ymb_hist import YMBHist

def H( x ):
    return ( sign( x ) + 1.0 ) / 2.0

class IDistrib( Interface ):
    pass

class YMBDistrib( YMBHist ):

    implements( IDistrib )

    n_int = Int( 30, auto_set=False, enter_set=True )

    x_array = Property( depends_on='slider.var_enum, slider.+changed_range, +modified' )
    @cached_property
    def _get_x_array( self ):
        h, b = self.normed_hist
        par_a = min( b )
        par_b = max( b )
        print par_a, par_b
        return linspace( par_a, par_b, self.n_int )

    pdf_array = Property( depends_on='slider.var_enum, slider.+changed_range, +modified' )
    @cached_property
    def _get_pdf_array( self ):

        # get the normed histogram implemented in the base class YMBHist
        h, b = self.normed_hist
        b_cen = ( b[1:] + b[:-1] ) / 2.
        mask = ( h != 0. )
        h = h[mask]
        b_cen = b_cen[mask]
        b = hstack( [min( b ), b_cen, max( b )] )
        h = hstack( [0, h, 0] )
        h = h / trapz( h, b )
        p = MFnLineArray( xdata=b, ydata=h )
        return p.get_values( self.x_array )

    data_mean = Property( Float, depends_on='slider, slider.var_enum, +modified' )
    @cached_property
    def _get_data_mean( self ):
        data = self.slider.stat_data
        return mean( data[data >= 0] )
    
    @on_trait_change( '+modified, slider.var_enum, slider.+changed_range, data.+changed_source, data.+changed_config' )
    def _redraw( self ):
        figure = self.figure
        axes = figure.axes[0]
        axes.clear()

        var_data = self.slider.stat_data
        
        axes.hist( var_data[var_data >= 0], bins=self.bins, normed=True )
        
        # @kelidas: plot the values of n_int values (it cannot agree with histogram) 
        axes.plot( self.x_array, self.pdf_array, color='red', linewidth=5 )
        axes.set_xlabel( self.slider.var_enum )
        
        if self.ylimit_on == True:
            axes.set_ylim( 0, self.ylimit )
        self.data_changed = True

if __name__ == '__main__':
    ymb_pd = YMBDistrib( varname='contact_fraction' )
    ymb_pd.configure_traits()

