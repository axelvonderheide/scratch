

from enthought.traits.api import HasTraits, Float, Property
from enthought.traits.ui.api import View, Item
from math import sqrt

class Foo( HasTraits ):
    
    parameter = Float( 10.0 )
    
    def _from_p_to_mom(self, p ):
        return p ** 2
    
    def _from_mom_to_p(self, m ):
        return sqrt(m)
    
    moment = Property( Float, depends_on = 'parameter' )
    def _get_moment(self):
        return self._from_p_to_mom( self.parameter )

    def _set_moment(self, value):
        self.parameter = self._from_mom_to_p(value)
        
    traits_view = View( Item('parameter'),
                        Item('moment' )
                        )
    
    
if __name__ == '__main__':
    f = Foo()
    f.configure_traits()