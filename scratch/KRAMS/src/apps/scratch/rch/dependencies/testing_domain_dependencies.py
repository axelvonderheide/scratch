
from enthought.traits.api import \
    HasTraits, Float, Int, Array, Property, cached_property, \
    Tuple, List, Str, on_trait_change

from enthought.traits.ui.api import \
    View, Item

class TestDomain( HasTraits ):
    
    a = Int(1)
    
    aa = Property( Int, depends_on = 'a' )
    @cached_property
    def _get_aa( self ):
        return self.a * self.a

    a_a = Property( Int, depends_on = 'a' )
    @cached_property
    def _get_a_a( self ):
        return self.a * self.a
        
    aa10 = Property( Int, depends_on = 'a' )
    @cached_property
    def _get_aa10(self):
        return self.a_a * 10
    
    traits_view = View(Item('a'),
                       Item('aa10', style = 'readonly' ),
                       resizable = True,
                       height = 0.5,width = 0.5 )

tp = TestDomain( a = 1 )

tp.configure_traits()
