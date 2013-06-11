
import sys
from enthought.traits.api import \
    HasTraits, Instance, Property, Int, cached_property
import gc

class OrdinaryClass(HasTraits):
    
    a = Int
        
    def __del__(self):
        print "deleting instance of ordinary class, a =", self.a

class ClassWithProperty(HasTraits):
    
    a = Int

    # The class gets deleted if this is not here
    a_plus_1 = Property( Int, depends_on ='a' )
    def _get_a_plus_1(self):
        return self.a + 1
    
    def __del__(self):
        # ?! THIS DOES NOT HAPPEN !?
        print "deleting instance of a class with dependent property"

def fn_scope():
    '''Function providing the scope for instances. 
    '''
    vanishing_object    = OrdinaryClass( a = 1 )
    print 'referrers to vanishing', gc.get_referrers( vanishing_object )
    persisting_object = ClassWithProperty( a = 1 )
    print 'referrers to persisting', gc.get_referrers( persisting_object )

# calling the function with local variables
#
fn_scope()  # if I put a loop here - the number of unreachables grows

# collecting the garbage
#
n_unreachable = gc.collect()

# there are 36 objects unreachable for the gc 
# per one instance of ClassWithProperty
#
print 'number of unreachable objects', n_unreachable

# the Heap browser shows one existing instance of ClassWithProperty
# with two referrers (one is probably the HeapBrowser itself.
# 
from enthought.developer.tools.heap_browser import HB_HeapBrowser
hb = HB_HeapBrowser()
hb.update()
hb.configure_traits()



