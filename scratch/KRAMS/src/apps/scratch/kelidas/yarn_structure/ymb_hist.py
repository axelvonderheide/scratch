#-------------------------------------------------------------------------------
#
# Copyright (c) 2009, IMB, RWTH Aachen.
# All rights reserved.
#
# This software is provided without warranty under the terms of the BSD
# license included in simvisage/LICENSE.txt and may be redistributed only
# under the conditions described in the aforementioned license.  The license
# is also available online at http://www.simvisage.com/licenses/BSD.txt
#
# Thanks for using Simvisage open source!
#
# Created on Dec 13, 2010 by: rch




from enthought.traits.api import HasTraits, Float, Property, cached_property, \
    Event, Array, Instance, File, Int, Directory, Button, Range, Enum, \
     on_trait_change, Bool, Trait, HasPrivateTraits, Constant, List, Tuple, \
     DelegatesTo
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
    sqrt, sum, all, zeros_like, zeros, ones, diff, where, unique, isnan, pi, invert
from numpy.random import random
from os.path import join
from promod.simdb import SimDB
from traits.editors.mpl_figure_editor import MPLFigureEditor
import matplotlib.pyplot as plt
import numpy.ma as ma
import os
import re
from enthought.traits.trait_types import DelegatesTo

from ymb_data import YMBData
from ymb_data import YMBSlider

class YMBHist( HasTraits ):
    #raw_data = Instance( YarnRawData, () )
    data = Instance( YMBData, () )

    slider = Instance( YMBSlider, () )

    figure = Instance( Figure, () )

    bins = Int( 20, auto_set=False, enter_set=True, modified=True )
    ylimit_on = Bool( False, modified=True )
    ylimit = Int( 100, auto_set=False, enter_set=True, modified=True )

    normed_hist = Property( depends_on='slider.var_enum, slider.+changed_range, +modified' )
    @cached_property
    def _get_normed_hist( self ):
        data = self.slider.stat_data
        h, b = histogram( data[data >= 0], bins=self.bins, normed=True )
        return h, b

    def _figure_default( self ):
        figure = Figure()
        figure.add_axes( [0.1, 0.1, 0.8, 0.8] )
        return figure

    data_changed = Event( True )
    @on_trait_change( '+modified, slider.var_enum, slider.+changed_range, data.+changed_source, data.+changed_config' )
    def _redraw( self ):
        figure = self.figure
        axes = figure.axes[0]
        axes.clear()
        var_data = self.slider.stat_data
        axes.hist( var_data[var_data >= 0], bins=self.bins )
        axes.set_xlabel( self.slider.var_enum )
        if self.ylimit_on == True:
            axes.set_ylim( 0, self.ylimit )
        self.data_changed = True

    view = View( HSplit( 
                Group( Item( 'slider', style='custom', show_label=False ),
                    label='yarn data',
                    id='yarn_hist.config',
                    dock='tab',
                ),
                    Group( Item( 'figure', style='custom',
                                  editor=MPLFigureEditor(),
                                  show_label=False ),
                    HGroup( 
                             Item( 'bins' ),
                             HGroup( 
                                    Item( 'ylimit_on', label='Y limit' ),
                                    Item( 'ylimit', enabled_when='ylimit_on == True',
                                          show_label=False ), )
                             ),
                    label='histogram',
                    id='yarn_hist.figure',
                    dock='tab',
                ),
                id='yarn_hist.split',
                dock='tab',
                ),
                id='yarn_structure_view',
                resizable=True,
                scrollable=True,
                dock='tab',
                width=0.8,
                height=0.4
                        )


if __name__ == '__main__':
    yarn = YMBHist()
    yarn.configure_traits()


