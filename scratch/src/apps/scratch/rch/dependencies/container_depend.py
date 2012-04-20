


from enthought.traits.api import \
     Array, Bool, Callable, Enum, Float, HasTraits, Interface, implements, \
     Instance, Int, Trait, Str, Enum, Callable, List, TraitDict, Any, \
     on_trait_change, Tuple, WeakRef, Delegate, Property, cached_property, Constant


from enthought.traits.ui.api \
    import TabularEditor, View, Item, Group, HGroup, VGroup, HSplit, VSplit 
from enthought.traits.ui.menu \
    import OKButton, CancelButton
from enthought.traits.ui.tabular_adapter \
    import TabularAdapter    
    

class Entry( HasTraits ):
    number = Float(20, affects_fn = True )
    number2 = Float(20 )
    

class EntryAdapter ( TabularAdapter ):
    '''This adapter specialization lists the attributes of RIDVariable
    to be displayed in columns of the TabularEditor'''
    
    columns = [ ( 'number',         'number' ),
               ( 'number2', 'number2' ) ]
                
    font                      = 'Courier 10'
    variable_alignment        = Constant( 'right' )
    
# -- Tabular Editor Definition -------------------------------------------------
# Edit the list of random variables
# 
entry_editor = TabularEditor(
    selected   = 'current_entry',
    adapter    = EntryAdapter(),
    operations = [ 'move' ],
    auto_update = True
)

class Container(HasTraits):
    
    cont = List( Entry )
    def _cont_default(self):
        return [Entry(),
                Entry(),
                Entry(),
                Entry(),
                ]
        
    current_entry = Instance( Entry )
    def _current_entry_default(self):
        return self.cont[0]
    
    @on_trait_change('cont.+affects_fn')
    def cont_updated(self):
        print 'cont_updated'
    
    traits_view = View( HSplit( 
                                   Item('cont', show_label = False, editor = entry_editor ),
                                   Item('current_entry', show_label = False, 
                                        style = 'custom', resizable = True),
                                   ),
                                   resizable = True,
                                   scrollable = True,
                                   height = 0.4,
                                   width = 0.6
                                   )

c = Container()
c.cont.append( Entry() )
c.configure_traits()
    