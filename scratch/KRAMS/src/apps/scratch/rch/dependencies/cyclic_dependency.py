

from enthought.traits.api import \
     Array, Bool, Callable, Enum, Float, HasTraits, Interface, implements, \
     Instance, Int, Trait, Str, Enum, Callable, List, TraitDict, Any, \
     on_trait_change, Tuple, WeakRef, Delegate, Property, cached_property

class A( HasTraits ):
    
    a = Float(0, auto_set = False, enter_set = True )
    
    b = Float(0, auto_set = False, enter_set = True )
    
    @on_trait_change('a')
    def change_b(self):
        print 'changing b due to a'
        self.on_trait_change( self.change_a, 'b', remove = True )
        self.b = self.a + 20
        self.on_trait_change( self.change_a, 'b' )
        
    @on_trait_change('b')
    def change_a(self): 
        print 'changing a due to b'
        self.on_trait_change( self.change_b, 'a', remove = True )
        self.a = self.b - 20
        self.on_trait_change( self.change_b, 'a' )
        
        
a = A()
a.configure_traits()