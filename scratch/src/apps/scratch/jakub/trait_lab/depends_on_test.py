'''
Created on May 18, 2009

@author: jakub
'''
# cached_prop.py -- Example of @cached_property decorator
from enthought.traits.api import HasPrivateTraits, List, Int,\
                                 Property, cached_property

class TestScores ( HasPrivateTraits ):

    scores  = List( Int )
    average = Property( depends_on = 'scores' )

    @cached_property
    def _get_average ( self ):
        s = self.scores
        return (float( reduce( lambda n1, n2: n1 + n2, s, 0 ) )
                 / len( s ))

    def __del__(self):
        print "deleting instance"
        

    ts = TestScores(scores =[1,2,3] )
    print "average ", ts.average
    
    import sys
    print 'this should be 2 but it is',sys.getrefcount(ts)
    
    print 'destructor never called'
    del ts