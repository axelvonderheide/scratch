

from enthought.traits.api import \
     Array, Bool, Callable, Enum, Float, HasTraits, \
     Instance, Int, Trait, Range, HasTraits, on_trait_change, Event, \
     implements, Dict, Property, cached_property, Str

from enthought.traits.ui.api import \
    View, Item

# Enthought library imports.
from enthought.pyface.api import ApplicationWindow, GUI
from enthought.pyface.action.api import Action, MenuManager, MenuBarManager
from enthought.pyface.action.api import StatusBarManager, ToolBarManager
from enthought.pyface.api import ApplicationWindow, SplitPanel


class A(HasTraits):
    print ' A called'
    val = Int(101)
    val2 = Int(102)
    

class B(HasTraits):
    
    a = Instance(A)
    def _a_default(self):
        print '_a_default called'
        return A()
    
    b = Float(2.)
                 
    c = Property(Float,depends_on = "a.+,b" )
#    c = Property(Float,depends_on = "a, b" )
    @cached_property
    def _get_c(self):
        print ' _get_c_called'
        return self.a.val * self.b * self.a.val2
    
    view = View( Item('a@'),
                 Item('b'),
                 Item('c'),
                 resizable = True )
    
print ' before A construction'
a = A()
print ' after A construction \n'

print ' before B construction'
b = B()
print ' after B construction'

print 'b.c', b.c
print 'b.c', b.c

b.configure_traits()