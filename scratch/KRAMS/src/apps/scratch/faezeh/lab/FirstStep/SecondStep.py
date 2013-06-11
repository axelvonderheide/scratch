
# Imports:
from enthought.traits.api \
    import HasTraits, HasPrivateTraits, Array, Tuple, Float, Button
    
from enthought.traits.ui.api \
    import View, HSplit, VSplit, Item, TupleEditor, Group
    
#from enthought.model.numeric_ui \
#    import Table, Plot

from math import cos

from enthought.traits.ui.menu import OKButton, CancelButton

from pylab import *


class AX2COSX(HasTraits):
    ''' Simple class emulating a function a * x**2 * cos(b*x+c)
    '''
    a = Float
    b = Float
    c = Float
    
    
    # Traits view definitions:  
    view = View(Group(Item( name = 'a'), 
                      Item( name = 'b'),
                      Item( name = 'c' ),
                      label = 'Numeric Rechnung', 
                      show_border = True))

    t = arange(0.0, 5.2, 0.2)

        
    def __call__(self, x, **kwargs ):
        return self.a*x**2 * cos(self.b*x + self.c)
    
if __name__ == '__main__':
    fx = AX2COSX( a = 1.0, b = 2.0, c = 3.0)
    v = fx(3.0, a = 20, some_b = 345 )
    print 'result', v
    fx.configure_traits()
