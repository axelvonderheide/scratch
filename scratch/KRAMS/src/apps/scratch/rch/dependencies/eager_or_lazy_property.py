
from enthought.traits.api import \
    Float, HasTraits, Str, Bool, List, Property, cached_property, Instance, \
    Constant, Button
    
class Simulator( HasTraits ):
    
    some_param = Float(1.0)
    show_some_param = Bool( True )
    
    input_param = Bool( True )
    
    response = Property( depends_on = 'input_param' )
    @cached_property
    def _get_response(self):
        print 'simulation running ...' 
        return 1

####### View #####
# -- Tabular Editor Definition -------------------------------------------------
from enthought.traits.ui.api \
    import View, Item, Group, HGroup, VGroup, HSplit, VSplit, Tabbed 

#
view = View(
            Item( 'show_some_param' ),
            Item( 'some_param', visible_when = 'object.show_some_param == True'),
            Item('input_param' ),
            title = 'Eager vs. Lazy',
            resizable = True,
            dock = 'tab',
            )

sim = Simulator()
sim.configure_traits( view = view ) # kind = 'modal' )
