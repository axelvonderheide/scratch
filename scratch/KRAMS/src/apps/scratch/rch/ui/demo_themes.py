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
# Created on Dec 17, 2009 by: rch

from enthought.traits.api import HasTraits, Str, Range, Float, Enum
from enthought.traits.ui.api import View, Group, Item, Label
from enthought.traits.ui.wx.themed_text_editor import \
    ThemedTextEditor

class Test ( HasTraits ):

    name   = Str
    age    = Range( 1, 100 )
    weight = Float
    gender = Enum( 'Male', 'Female' )

    view = View(
        Group(
            Group(
                Label( 'A Themed Label', '@GF6' ),
                Item( 'name' ),
                Item( 'age' ),
                Item( 'weight', editor=ThemedTextEditor()),
                Item( 'gender' ),
                group_theme = '@GD0'
            ),
            group_theme = '@G',
            item_theme  = '@B0B',
            label_theme = '@BEA'
        ),
        title   = 'Themed Traits UI',
    )

Test().configure_traits()
