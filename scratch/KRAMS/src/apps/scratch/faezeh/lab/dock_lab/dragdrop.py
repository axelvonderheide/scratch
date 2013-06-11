from enthought.traits.api \
import HasTraits, Any, Str

from enthought.traits.ui.api import View, Item
    
class DroppedObjectValue ( HasTraits ):
# Holds the most recently dropped object:
    object = Any( droppable='Drop any object here to see '
                      'its value' )
# Contains string value of the most recently dropped
# object:
    value=Str('Drop an object on me to display its value')
    traits_view = View(
                       Item( 'value',
                             style
                             = 'readonly',
                             show_label = False
                             )
                       )
# Handle a new object being dropped:
def _object_changed ( self, object ):
    self.value = str( object )
# Create an importable application object:
dropped_object_value = DroppedObjectValue()

if __name__ == '__main__':
    dropped_object_value.configure_traits()