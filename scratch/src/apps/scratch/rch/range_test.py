

from enthought.traits.api import Range, HasTraits, Trait, Button, Float
from enthought.traits.ui.api import View, Item, RangeEditor

class RangeContainer(HasTraits):
    
    first_range = Range( 0., 2000., 10. )
    alpha = Trait(1.0, Range(0.0, 1.0))


    start_time = Float( 0. )
    stop_time = Float( 1. )
    delta_time = Float( 0.01 )
    current_time = Float( 0. )

    run_loop = Button()
    def _run_loop_fired(self):
        self.current_time += self.delta_time
    
    traits_view = View( Item( 'start_time' ),
                        Item( 'stop_time' ),
                        Item( 'delta_time' ), 
                        Item( 'current_time', style = 'custom',
                              editor = RangeEditor( low_name = 'start_time',
                                                    high_name = 'stop_time',
                                                    mode = 'slider' ) ),
                        Item('run_loop'),
                        Item('alpha'),
                        resizable = True )
    

rc = RangeContainer()
rc.first_range = 12
rc.alpha = 0.4
rc.configure_traits()
    
