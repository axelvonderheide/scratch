
from enthought.traits.api import \
    HasTraits, Float, Int, Array, Property, cached_property, \
    Tuple, List, Str, on_trait_change

from enthought.traits.ui.api import \
    View, Item

class TestProperty( HasTraits ):
    
    call_history = List( Str )
    call_counter = 0

    # the source attribute - starting option everything depends on it
    myattrib = Str('default', auto_set = False, enter_set = True )
    def _myattrib_default(self):
        self.call_history.append('_myattrib_default')
        return 'default[' + `self.call_counter` + ']'

    # myattrib -> prop1
    prop1 = Property( Str, depends_on = 'myattrib' )
    @cached_property
    def _get_prop1(self):
        self.call_history.append('_get_prop1')
        self.call_counter += 1
        s = 'prop1[' + str(self.call_counter) + ']'
        return s + '( ' + self.myattrib + ' )'

    # myattrib -> prop2
    prop2 = Property( Str, depends_on = 'myattrib' )
    @cached_property
    def _get_prop2(self):
        self.call_history.append('_get_prop2')
        self.call_counter += 1
        s = 'prop2[' + str(self.call_counter) + ']'
        return s + '( ' + self.myattrib + ' )'
        
    # myattrib -> prop1 -> prop2
    prop3 = Property( Str, depends_on = 'myattrib' )
    @cached_property
    def _get_prop3(self):
        self.call_history.append('_get_prop3')
        self.call_counter += 1
        s = 'prop3[' + str(self.call_counter ) + ']'
        return s + '( ' + self.prop1 + ' )'

    def access_property(self):
        return  'access: ' + self.prop3

    traits_view = View(Item('myattrib'),
                            'prop1','prop2','prop3','call_history@',
                       resizable = True,
                       height = 0.5,width = 0.5 )

tp = TestProperty( myattrib = 'init' )

print 'ACCESSING prop2'
print '>>>', tp.prop2

print 'CHANGING myattrib'
tp.myattrib = 'change1'

print 'ACCESSING prop2'
print '>>>',tp.prop2

print 'CHANGING myattrib'
tp.myattrib = 'change2'
print '>>>',tp.access_property()

tp.configure_traits()

