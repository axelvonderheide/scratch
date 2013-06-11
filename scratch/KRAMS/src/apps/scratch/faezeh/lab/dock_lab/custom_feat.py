from enthought.traits.api \
    import HasTraits, Int, List
from enthought.traits.ui.api import  View, Item
from enthought.pyface.image_resource \
    import ImageResource
#from enthought.pyface.dock.features.api \
#    import CustomFeature
from enthought.traits.api \
import *


class CustomFeature ( HasPrivateTraits ):
# The current image to display on the DockControl tab:
    tab_image = Instance( ImageResource )
    # The current (optional) image to display on the
    # DockControl drag bar:
    bar_image = Instance( ImageResource )
    # The tooltip to display when the mouse is hovering
    # over the image:
    tooltip = Str
    # Is the feature currently enabled?
    enabled = true
    # Name of the method to invoke on a left click:
    left_click = Str
    # Name of the method to invoke on a right click:
    right_click = Str
    # Name of the method to invoke when the user starts to
    # drag:
    drag = Strcustom_feature = True
    # Name of the method to invoke when the user starts to
    # ctrl-drag:
    control_drag = Str
    # Name of the method to invoke when the user starts to
    # shift-drag:
    shift_drag = Str
    # Name of the method to invoke when the user starts to
    # alt-drag:
    alt_drag = Str
    # Name of the method to invoke when the user drops an
    # object:
    drop = Str
    # Name of the method to invoke to see if the user can
    # drop an object:
    can_drop = Str
    
class Counter ( HasTraits ):
# The current count value:
    count = Int
    
    IncrementFeature = CustomFeature(
                                 tab_image = ImageResource( 'close' ),
                                 tooltip = 'Click to increment count',
                                 left_click = 'increment')
    DecrementFeature = CustomFeature(
                                 tab_image = ImageResource( 'decrement' ), 
                                 tooltip = 'Click to decrement count',
                                 left_click = 'decrement')

    # The feature definitions:
    features = List( [IncrementFeature, DecrementFeature],
                     custom_feature = True )
    # The application object view:
    traits_view = View(Item( 'count', style = 'readonly'),
                       #Item( 'features', style = 'readonly'),
                       title="using custom feature",
                       resizable=True,
                       dock = 'tab',
                       id = 'myfeature',
                       width=400, height=300)
    # Methods to handle 'feature' left click actions:
    def increment ( self ): 
        self.count += 1
    def decrement ( self ): 
        self.count -= 1
# Create an importable instance:
#counter = Counter()

if __name__ == '__main__':
    Counter().configure_traits()
