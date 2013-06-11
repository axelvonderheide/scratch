'''
Created on May 3, 2010

@author: alexander
'''
from enthought.traits.api import \
    PrototypedFrom, Float, HasTraits, Instance, Str, DelegatesTo



class Parent ( HasTraits ):

    first_name = Str
    last_name  = Str

    def _last_name_changed(self, new):
        print "Parent's last name changed to %s." % new

class Child ( HasTraits ):

    father = Instance( Parent )
    first_name = Str
#    last_name  = DelegatesTo( 'father' )
#    last_name  = DelegatesTo( 'father', listenable = False )
    last_name  = PrototypedFrom( 'father' )

    def _last_name_changed(self, new):
        print "Child's last name changed to %s." % new

dad = Parent( first_name='William', last_name='Chase' )
son = Child( first_name='John', father=dad )
print 'dad.last_name', dad.last_name
print 'son.last_name', son.last_name

dad.last_name='Jones'

print 'dad.last_name', dad.last_name
print 'son.last_name', son.last_name

son.last_name='Thomas'

print 'dad.last_name', dad.last_name
print 'son.last_name', son.last_name

dad.last_name='Riley'

print 'dad.last_name', dad.last_name
print 'son.last_name', son.last_name

del son.last_name

print 'dad.last_name', dad.last_name
print 'son.last_name', son.last_name

dad.last_name='Simmons'

print 'dad.last_name', dad.last_name
print 'son.last_name', son.last_name
