'''
Created on Nov 20, 2010

@author: kelidas
'''
#  Copyright (c) 2007, Enthought, Inc.
#  License: BSD Style.
#-- Imports --------------------------------------------------------------------

import logging, sys
logging.basicConfig( stream=sys.stderr )

from random \
    import choice
    
from enthought.traits.api \
    import HasPrivateTraits, Str, Enum, Range, List, Button, Instance, \
           Property, cached_property, on_trait_change
    
from enthought.traits.ui.api \
    import View, VGroup, HGroup, Item, ListEditor, spring

#-- The Hotel class ------------------------------------------------------------

class Hotel ( HasPrivateTraits ):
    
    # The season of the year:
    season = Enum( 'Winter', 'Spring', 'Summer', 'Fall' )
    
    # The current cost of heating fuel (in dollars/gallon):
    fuel_cost = Range( 2.00, 10.00, 4.00 )
    
    # The current minimum temparature allowed by the hotel:
    min_temperature = Property( depends_on='season, fuel_cost' )
    
    # The guests currently staying at the hotel:
    guests = List # ( Instance( 'Guest' ) )
    
    # Add a new guest to the hotel:
    add_guest = Button( 'Add Guest' )
    
    # The view of the hotel:
    view = View( 
        VGroup( 
            HGroup( 
                Item( 'season' ), '20',
                Item( 'fuel_cost', width=300 ),
                spring,
                Item( 'add_guest', show_label=False ),
                show_border=True,
                label='Hotel Information'
            ),
            VGroup( 
                Item( 'guests',
                      style='custom',
                      editor=ListEditor( use_notebook=True,
                                           deletable=True,
                                           dock_style='tab',
                                           page_name='.name' )
                ),
                show_labels=False,
                show_border=True,
                label='Guests'
            )
        ),
        title='The Belmont Hotel Dashboard',
        width=0.6,
        height=0.2,
        resizable=True
    )
    
    # Property implementations:
    @cached_property
    def _get_min_temperature ( self ):
        return ( { 'Winter': 32,
                  'Spring': 40,
                  'Summer': 45,
                  'Fall':   40 }[ self.season ] + 
                  min( int( 60.00 / self.fuel_cost ), 15 ) )
        
    # Event handlers:
    @on_trait_change( 'guests[]' )
    def _guests_modified ( self, removed, added ):
        for guest in added:
            guest.hotel = self
            
    def _add_guest_changed ( self ):
        self.guests.append( Guest() )

#-- The Guest class ------------------------------------------------------------

class Guest ( HasPrivateTraits ):
    
    # The name of the guest:
    name = Str
    
    # The hotel the guest is staying at:
    hotel = Instance( Hotel )
    
    # The room plan the guest has chosen:
    plan = Enum( 'Flop house', 'Cheap', 'Cozy', 'Deluxe' )
    
    # The maximum temperature allowed by the guest's plan:
    max_temperature = Property( depends_on='plan' )
    
    # The current room temperature as set by the guest:
    temperature = Range( 'hotel.min_temperature', 'max_temperature' )
    
    # The view of the guest:
    view = View( 
        Item( 'plan' ),
        Item( 'temperature' )
    )
    
    # Property implementations:
    @cached_property
    def _get_max_temperature ( self ):
        return { 'Flop house': 62,
                 'Cheap':      66,
                 'Cozy':       75,
                 'Deluxe':     85 }[ self.plan ]
                 
    # Default values:
    def _name_default ( self ):
        return choice( 
            [ 'Leah', 'Vibha', 'Janet', 'Jody', 'Dave', 'Evan', 'Ilan', 'Gael',
              'Peter', 'Robert', 'Judah', 'Eric', 'Travis', 'Mike', 'Bryce',
              'Chris' ] )
    
#-- Create the demo ------------------------------------------------------------

# Create the demo object:
demo = Hotel( guests=[ Guest() for i in range( 5 ) ] )

# Run the demo (if invoked from the command line):
if __name__ == '__main__':
    logging.info( 'Start!' )
    demo.configure_traits()
