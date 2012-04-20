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
# Created on Jan 19, 2010 by: rch

from enthought.traits.api import HasTraits, List, Int
from enthought.traits.ui.api import View, Item, CheckListEditor, EnumEditor

class foo(HasTraits):
    thelist = List()
    def _thelist_default(self):
        print 'here'
        return [self.items[1]]

    items = List()

    traits_view = View(Item(name = 'thelist',
                            editor = EnumEditor( name = 'items' ),
                            style = 'simple', show_label = False
                            )
                      )

f = foo(items=[ (1, 'one'), (2, 'two')])
print f.thelist  # before edit
f.configure_traits()
print f.thelist  # after edit
