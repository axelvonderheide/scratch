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
# Created on Dec 1, 2009 by: rch

from enthought.traits.api import *
from enthought.traits.ui.api import *
from enthought.traits.ui.table_column import *

class Person( HasTraits ):
    pass

class Parent ( Person ):
    
    first_name = Str
    last_name  = Str
    
class Child ( Person ):
    
    mother = Instance( Parent )
    father = Instance( Parent )
    
    first_name = Str
    last_name  = Delegate( 'father' )
    
class ChildModelView ( ModelView ):
    
    # Define the 'family' ModelView property that maps the child and its 
    # parents into a list of objects that can be viewed as a table:
    family = Property( List )
    
    selected = Instance( Person )
    
    # Define a view showing the family as a table:
    view = View(
        Item( 'family', 
              show_label = False,
              editor = TableEditor(
                  auto_add = False,
                  selection_mode = 'row',
                  row_factory = Parent,
                  selected = 'selected',
                  editable = True,
                  show_toolbar = True,
                  deletable = True, 
                  columns = [ ObjectColumn( name = 'first_name' ),
                              ObjectColumn( name = 'last_name' ) ] ) ),
        Item( 'selected@', show_label = False ),
        resizable = True,
        height = 0.5,
        width = 0.5,
        buttons = ['OK','Cancel' ]
    )
       
    # Implementation of the 'family' property:
    def _get_family ( self ):
        return [ self.model.father, self.model.mother, self.model ]
    
# Create a sample family:
mom = Parent( first_name = 'Julia', last_name = 'Wilson' )
dad = Parent( first_name = 'William', last_name = 'Chase' )        
son = Child( mother = mom, father = dad, first_name = 'John' )

# Create the controller for the model:
demo = ChildModelView( model = son )
demo.configure_traits()
