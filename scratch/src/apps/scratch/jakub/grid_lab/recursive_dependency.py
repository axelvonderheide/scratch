'''
Created on May 4, 2009

@author: jakub
'''

from enthought.traits.api import HasTraits, This, Bool, Int, Str, Property, cached_property

class SubDomain( HasTraits ):
    
    name = Str
    previous = This
    state_changed = Bool(False)
    
    size = Int( 10 )
    
    @on_trait_change
    def size_changed(self):
        state_changed = True
    
    offset = Property( depends_on = 'previous.state_changed' )
    @cached_property
    def _get_offset(self):
        print '*** recalculating offset', self.name
        if self.previous == None:
            return 0
        else:
            return self.previous.offset + self.previous.size
    
print '>>> d1 = SubDomain()'
d1 = SubDomain( name = 'd1' )
print '<<< d1 = SubDomain()'
print '>>> d2 = SubDomain( previous = d1 )'
d2 = SubDomain( name = 'd2', previous = d1 )
print '<<< d2 = SubDomain( previous = d1 )'
print '>>> d3 = SubDomain( previous = d2 )'
d3 = SubDomain( name = 'd3', previous = d2 )
print '<<< d3 = SubDomain( previous = d2 )'
print '>>> d4 = SubDomain( previous = d3 )'
d4 = SubDomain( name = 'd4', previous = d3 )
print '<<< d4 = SubDomain( previous = d3 )'

print '>>> print d3.offset'
print d3.offset
print '<<< print d3.offset'

print '>>> d2.size = 5'
d2.size = 5
print '<<< d2.size = 5'

print '>>> print d3.offset'
print d3.offset
print '<<< print d3.offset'

print '>>> d1.size = 5'
d1.size = 5
print '<<< d1.size = 5 (eager call performed)'

print '>>> print d4.offset'
print d4.offset
print '<<< print d4.offset'
