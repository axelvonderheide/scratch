'''
Created on May 14, 2009

@author: jakub
'''
from enthought.traits.api import Array, Bool, Enum, Float, HasTraits, \
                                 Instance, Int, Trait, Str, Enum, \
                                 Callable, List, TraitDict, Any, Range, \
                                 Delegate, Event, on_trait_change, Button, Property, \
                                 cached_property, property_depends_on
from enthought.traits.ui.api import \
    Item, View, HGroup, ListEditor, VGroup, \
    HSplit, Group, Handler, VSplit, RangeEditor, spring
    
class TLine(HasTraits):
    '''
    Time line for the control parameter.
    
    This class sets the time-range of the computation - the start and stop time.
    val is the value of the current time.
    
    TODO - the info page including the number of load steps 
    and estimated computation time.
    
    TODO - the slide bar is not read-only. How to include a real progress bar?
    '''
    min = Float( 0.0 )
    max = Float( 1.0 )
    step = Float( 0.1 )
    val = Float( 0.0 )    
    mode = Enum( 'auto', 'slider', 'xslider', 'spinner', 'enum', 'text', 'logslider' )
    #time = Property
    traits_view = View( HGroup( Item('min'), 
                                spring, 
                                Item('step', label = 'step size'), 
                                spring, 
                                Item('max'),
                                spring,
                                Item('mode')), 
                        Item('time', editor = RangeEditor( low_name = 'min',
                                                           high_name = 'max',
                                                           format = '(%s)',
                                                           auto_set = False,
                                                           enter_set = False,
                                                           ),
                                    show_label = False
                                                           ),
                        resizable = True
                        )
#mode = Enum( 'auto', 'slider', 'xslider', 'spinner', 'enum', 'text', 'logslider' )
if __name__ == '__main__':
    tl = TLine()
    tl.configure_traits()