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
# Created on Mar 3, 2011 by: rch

from enthought.traits.api import \
    HasTraits, Enum, Instance
from enthought.traits.ui.api import \
    View, Item, Group, VGroup, HGroup, Tabbed

from stats.spirrid import RFModelView

class RFBrowser( HasTraits ):

    rf = Enum( values = 'rf_values' )
    def _rf_default( self ):
        return self.rf_values[0]
    def _rf_changed( self ):
        self.rf_model_view = RFModelView( model = self.rf )

    rf_model_view = Instance( RFModelView )
    def _rf_model_view_default( self ):
        return RFModelView( model = self.rf )

    traits_view = View( 
                  VGroup( 
                   Tabbed( 
                   VGroup( 
                         Item( 'rf', show_label = False ),
                         Item( 'rf_model_view@', show_label = False, resizable = True ),
                         label = 'Deterministic model',
                         id = 'spirrid.tview.model',
                         ),
                    scrollable = True,
                    id = 'spirrid.tview.tabs',
                    dock = 'tab',
            ),
            ),
            title = 'RF-Browser',
            id = 'spirrid.viewmodel',
            dock = 'tab',
            resizable = True,
            height = 1.0, width = 1.0,
            buttons = ['OK', 'Cancel']
            )

if __name__ == '__main__':

    from stats.spirrid.rf_filament import \
        Filament

    # in the crack bridge module

    from quaducom.crackbridge.double_sided_yarn import \
        DoublePullout
    from quaducom.crackbridge.crack_bridge import \
        StressInFiberWithConstantFriction
    from quaducom.crackbridge.yarn_symmetrical import \
        DoublePulloutSym

    # in the pullout module

    from quaducom.pullout.const_frict_free_length_finite_fiber import \
        ConstantFrictionFreeLengthFiniteFiber
    from quaducom.pullout.constant_friction_finite_fiber import \
        ConstantFrictionAndFreeLength, ConstantFrictionFiniteFiber
    from quaducom.pullout.double_sidded_yarn import \
        DoublePullout as DP

    # Currently available response functions.
    #

    rf_list = [
               Filament(),
               StressInFiberWithConstantFriction(),
               DoublePullout(),
               DoublePulloutSym(),
               ConstantFrictionFreeLengthFiniteFiber(),
               ConstantFrictionAndFreeLength(),
               ConstantFrictionFiniteFiber(),
               DP()
               ]

    rfb = RFBrowser( rf_values = rf_list )
    rfb.configure_traits()
