
from enthought.traits.api import HasTraits, This, Int, Property, cached_property

class SubDomain( HasTraits ):
    
    previous = This
    
    size = Int( 10, changed_state = True )
    
    offset = Property( depends_on = 'previous.+changed_state', changed_state = True )
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
print d4.offset # show the offset of the last domain
print '<<< print d3.offset'

#print '>>> d2.size = 5'
#d2.size = 5 # reduce the size of the second domain
#print '<<< d2.size = 5'
#
#print '>>> print d3.offset'
#print d3.offset # look at the offset of the third domain
#print '<<< print d3.offset'

print '>>> d1.size = 5'
d1.size = 5 # reduce the size of the first subdomain
print '<<< d1.size = 5 (eager calls in previous domains performed)'

print '>>> print d3.offset'
print d4.offset # look at the offset of the third domain
print '<<< pint d3.offset'
