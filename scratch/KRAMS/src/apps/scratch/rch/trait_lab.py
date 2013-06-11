

from enthought.traits.api import HasTraits, Int, cached_property, Property, \
    Interface, implements, Instance, Delegate

class IB( Interface ):
    pass
    #my_int = Int
    
class B( HasTraits ):
    implements(IB)
    my_int = Int(10)

class A( HasTraits ):
    
    some_b = Instance( IB )
    
    my_int = Delegate('some_b')
    
    pint = Property( Int )
    _pint = Int( None )
    def _set_pint(self, value):
        self._pint = value
    def _get_pint(self):
        if self._pint == None:
            return 200
        else:
            return self._pint
    
a = A( some_b = B() )
print a.pint
print a.my_int
