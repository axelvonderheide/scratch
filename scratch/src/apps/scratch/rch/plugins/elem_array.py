
#-- Imports --------------------------------------------------------------------

from os.path \
    import join, dirname
    
from numpy \
    import sqrt
    
from numpy.random \
    import random

from enthought.traits.api \
    import HasTraits, Property, Array, Any, Event, \
    on_trait_change, Int
    
from enthought.traits.ui.api \
    import View, Item, TabularEditor
    
from enthought.traits.ui.menu \
    import NoButtons, CancelButton
    
from enthought.traits.ui.tabular_adapter \
    import TabularAdapter

from enthought.pyface.image_resource \
    import ImageResource
    
#-- Constants ------------------------------------------------------------------

# Necessary because of the dynamic way in which the demos are loaded:
import enthought.traits.ui.api

#-- Tabular Adapter Definition -------------------------------------------------

class ArrayAdapter ( TabularAdapter ):

    columns = Property
    def _get_columns(self):
        n_columns = getattr( self.object, self.name ).shape[1]
        cols = [ (str(i), i ) for i in range(3) ]
        return [ ('element', 'index') ] + cols
        
    #columns = [ ( 'i', 'index' ), ( 'x', 0 ), ( 'y', 1 ),  ( 'z', 2 ) ]
                
    font        = 'Courier 10'
    alignment   = 'right'
    format      = '%d'
    index_text  = Property
#    index_image = Property
    
    def _get_index_text ( self ):
        return str( self.row )
        
    def x_get_index_image ( self ):
        x, y, z = self.item
        if sqrt( (x - 0.5) ** 2 + (y - 0.5) ** 2 + (z - 0.5) ** 2 ) <= 0.25:
            return 'red_flag'
        return None

#-- Tabular Editor Definition --------------------------------------------------

tabular_editor = TabularEditor(
    adapter = ArrayAdapter(),
    selected_row = 'current_row',
#    operations = [ 'move' ],
#    auto_update = True
    
#    images  = [ ImageResource( 'red_flag', search_path = search_path ) ]
)

#-- ShowArray Class Definition -------------------------------------------------

class ShowArray ( HasTraits ):

    data = Array
    
    current_row = Int(-1)
    
    @on_trait_change('current_row')
    def _display_current_row(self):
        print 'row value'
        print self.current_row
    
    view = View(
        Item( 'data', editor = tabular_editor, show_label = False, style = 'readonly' ),
        title     = 'Array Viewer',
        width     = 0.6,
        height    = 0.5,
        resizable = True,
        buttons   = [CancelButton]
    )
    
# Run the demo (if invoked from the command line):
if __name__ == '__main__':

    from enthought.traits.api import Button
    class Container(HasTraits):
        show_array = Button
        def _show_array_fired(self):
            # Create the demo:
            demo = ShowArray( data = random( ( 100000, 3 ) ) )
            demo.configure_traits()
        view = View('show_array')
    
    c = Container()
    c.configure_traits()
    