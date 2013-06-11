
'''
This example demonstrates the dependency management.
'''


from enthought.traits.api import \
    HasTraits, Float, Int, Array, Property, cached_property, \
    Tuple, List, Str, on_trait_change

from enthought.traits.ui.api import \
    View, Item

class TestProperty( HasTraits ):
    
    call_history = List( Str )
    call_counter = 0

    # the source attribute - starting option everything depends on it
    myattrib = Str('default' )
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

    # myattrib -> prop1 -> prop2
    prop2 = Property( Str, depends_on = 'myattrib' )
    @cached_property
    def _get_prop2(self):
        self.call_history.append('_get_prop2')
        self.call_counter += 1
        s = 'prop2[' + str(self.call_counter ) + ']'
        return s + '( ' + self.prop1 + ' )'

    # invoked when prop2 has changed
    @on_trait_change('myattrib')
    def launch_refresh(self):
        self.call_history.append('launch_refresh')
        print 'called using on_trait_change'
        print self.myattrib
        print self.prop2
        
#    traits_view = View(Item('myattrib'),
#                            'prop1','prop2','prop3','call_history@',
#                       resizable = True,
#                       height = 0.5,width = 0.5 )

t_p = TestProperty( myattrib = 'init' )

print 'ACCESSING prop2 with init'
print '>>>', t_p.prop2

print 'CHANGING myattrib to change1'
t_p.myattrib = 'change1'

print 'ACCESSING prop2'
print '>>>', t_p.prop2
#
for call in t_p.call_history:
    print call

#tp.configure_traits()

