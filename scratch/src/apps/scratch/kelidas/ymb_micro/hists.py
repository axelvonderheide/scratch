'''
Created on May 18, 2011

@author: kelidas
'''

from quaducom.ymb_micro.ymb_data import YMBCutData, YMBSegmentData, YMBSlider, YMBSource, yarn_list, var_dict
from numpy import linspace, array, searchsorted

from enthought.traits.api import HasTraits, Property, cached_property, Event, \
    Instance, Int, on_trait_change, Bool, Str, Tuple, List, Array, Enum, Trait, Range
from enthought.traits.ui.api import Item, View, Group, HGroup, HSplit
from traits.editors.mpl_figure_editor import MPLFigureEditor
from matplotlib.figure import Figure

cf_limit_arr = linspace( 0, 1, 5 )

class Hists_data( HasTraits ):
    cf_limit_enum = Enum( list( cf_limit_arr ) , modified = True )
    yarn_type = Enum( yarn_list, modified = True )
    var_enum = Trait( 'slack', var_dict, modified = True )

    source = Property( Instance( YMBSource ), depends_on = 'yarn_type' )
    @cached_property
    def _get_source( self ):
        return YMBSource( yarn_type = self.yarn_type )

    hists_data = Property( Array, depends_on = 'source, var_enum, yarn_type' )
    @cached_property
    def _get_hists_data( self ):
        data_list = []
        for cf_limit in cf_limit_arr:
            data = YMBCutData( source = self.source, cf_limit = cf_limit )
            slider = YMBSlider( var_enum = self.var_enum, data = data )
            data_list.append( slider.stat_data )
        return data_list

    trats_view = View( 'yarn_type',
                       Item( 'cf_limit_enum', style = 'custom' ),
                       'var_enum' )

class Hists( HasTraits ):
    data = Instance( Hists_data )

    bins = Int( 30, auto_set = False, enter_set = True, modified = True )

    figure = Instance( Figure )
    def _figure_default( self ):
        figure = Figure()
        figure.add_axes( [0.15, 0.15, 0.75, 0.75] )
        return figure

    data_changed = Event( True )
    @on_trait_change( 'data.+modified, +modified' )
    def _redraw( self ):
        figure = self.figure
        axes = figure.axes[0]
        axes.clear()

        axes.hist( self.data.hists_data, bins = self.bins, histtype = 'step', linewidth = 1, \
                        edgecolor = '0.5' )
        axes.set_title( self.data.var_enum )
        axes.hist( self.data.hists_data[searchsorted( cf_limit_arr, self.data.cf_limit_enum )],
                   bins = self.bins, histtype = 'step', linewidth = 2, \
                    edgecolor = 'black', label = 'cf_limit = %s' % self.data.cf_limit_enum )
        axes.legend()
        # redefine xticks labels
        #inter = axes.xaxis.get_view_interval()
        #axes.set_xticks( linspace( inter[0], inter[1], 5 ) )
        #axes.set_xlabel( self.slider.var_enum )
        #axes.set_ylabel( 'frequency' )#, fontsize = 16
        self.data_changed = True

    traits_view = View( 
                       Item( 'data@', show_label = False ),
                       'bins',
                       Group( Item( 'figure', style = 'custom',
                              editor = MPLFigureEditor(),
                              show_label = False )
                              , id = 'figure.view' ),
                    resizable = True,
                    title = 'Histograms',
                     )

if __name__ == '__main__':
    data = Hists_data()
    hists = Hists( data = data )
    hists.configure_traits()



#exit()
#
#
#import matplotlib.pyplot as plt
#
#
#cf_limit_arr = linspace( 0, 1, 11 )
#
#yarn_type = 'MAG'
#
#source = YMBSource( yarn_type = yarn_type )
#
#for cf_limit in cf_limit_arr:
#    data = YMBCutData( source = source, cf_limit = cf_limit )
#    plt.figure( 0 )
#    slider = YMBSlider( var_enum = 'slack', data = data )
#    plt.title( slider.var_enum )
#    plt.hist( slider.stat_data, bins = 30, histtype = 'step', linewidth = 2, label = 'cf = %s' % cf_limit )
#    plt.legend( loc = 0 )
#
#    plt.figure( 1 )
#    slider = YMBSlider( var_enum = 'bond free length', data = data )
#    plt.title( slider.var_enum )
#    plt.hist( slider.stat_data, bins = 30, histtype = 'step', linewidth = 2, label = 'cf = %s' % cf_limit )
#    plt.legend( loc = 9 )
#
#plt.show()


