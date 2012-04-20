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

import wxversion
wxversion.select( '2.8' )

from enthought.traits.api import \
    HasTraits, Instance, Property, cached_property

from enthought.traits.ui.api import \
    View, Item, Group, HSplit, VSplit, HGroup, VGroup, Tabbed

from ymb_data import YMBData
from ymb_data import YMBSlider
from ymb_hist import YMBHist
from ymb_correl import YMBCutCorrelView, YMBCutCorrel
from ymb_view3d import YMBView3D
from ymb_view2d import YMBView2D

class YMBMicro( HasTraits ):
    '''
    Yarn-Matrix-Bond Microstructure analysis
    '''

    yarn_data = Instance( YMBData )
    def _yarn_data_default( self ):
        return YMBData()

    view_3d = Property( Instance( YMBView3D ), depends_on = 'yarn_data.input_change' )
    @cached_property
    def _get_view_3d( self ):
        return YMBView3D( data = self.yarn_data )

    cut_view = Property( Instance( YMBView2D ), depends_on = 'yarn_data.input_change' )
    @cached_property
    def _get_cut_view( self ):
        return YMBView2D( data = self.yarn_data )

    slider = Property( Instance( YMBSlider ), depends_on = 'yarn_data.input_change' )
    @cached_property
    def _get_slider( self ):
        return YMBSlider( data = self.yarn_data )

    histogram = Property( Instance( YMBHist ), depends_on = 'yarn_data.input_change' )
    @cached_property
    def _get_histogram( self ):
        return YMBHist( slider = self.slider )

    correlation_data = Property( Instance( YMBCutCorrel ), depends_on = 'yarn_data.input_change' )
    @cached_property
    def _get_correlation_data( self ):
        return YMBCutCorrel( yarn_data = self.yarn_data )

    correlation = Property( Instance( YMBCutCorrelView ), depends_on = 'yarn_data.input_change' )
    @cached_property
    def _get_correlation( self ):
        return YMBCutCorrelView( yarn_data = self.yarn_data,
                                  correl_data = self.correlation_data )

    traits_view = View( VSplit( 
                               Group( 
                                     Item( 'yarn_data@', show_label = False ),
                                     label = 'Data source'
                                     ),
                               Tabbed( 
                                   HSplit( 
                                       Group( 
                                             Item( 'view_3d@', show_label = False ),
                                             label = '3D view',
                                             id = 'qdc.ymb.view3d',
                                             dock = 'tab',
                                             ),
                                       Group( 
                                             Item( 'cut_view@', show_label = False ),
                                             label = '2D view',
                                             id = 'qdc.ymb.view2d',
                                             dock = 'tab',
                                             ),
                                    dock = 'tab',
                                    id = 'qdc.ymb.data_verification',
                                    label = 'Data verification',
                                    ),
                                   HSplit( 
                                       Group( 
                                             Item( 'histogram@', show_label = False ),
                                             label = 'Histogram',
                                             id = 'qdc.ymb.histogram',
                                             dock = 'tab',
                                             ),
                                       Group( 
                                             Item( 'correlation@', show_label = False ),
                                             label = 'Correlation',
                                             id = 'qdc.ymb.correlation',
                                             dock = 'tab',
                                             ),
                                    dock = 'tab',
                                    id = 'qdc.ymb.data_statistics',
                                    label = 'Statistics',
                                    ),
                                 dock = 'tab',
                                 id = 'gdc.ymb.tasks',
                                 ),
                               dock = 'tab',
                               id = 'qdc.ymb.split',
                               ),
                        id = 'qdc.ymb',
                        resizable = True,
                        scrollable = True,
                        dock = 'tab',
                        width = 0.8,
                        height = 0.4
                        )

if __name__ == '__main__':
    ymbmicro = YMBMicro()
    ymbmicro.configure_traits()
