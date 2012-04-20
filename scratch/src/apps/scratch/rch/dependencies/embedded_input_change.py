'''
Created on Apr 30, 2010

@author: alexander
'''

from enthought.traits.api import \
    HasTraits, Instance, Int, on_trait_change, Event

import os

class A( HasTraits ):
    
    value = Int( 0, input = True )

    input_change = Event()
    @on_trait_change('+input')
    def _set_input_change(self):
        print 'A: input change raised'
        self.input_change = True
    
class B( HasTraits ):

    a = Instance( A )
    def _a_default(self):
        return A()
    
    input_change = Event()
    @on_trait_change('+input,a.input_change')
    def _set_input_change(self):
        print 'B: input change propagated'
        self.input_change = True

class C( HasTraits ):
    
    b = Instance( B )
    def _b_default(self):
        return B()
    
    input_change = Event()
    @on_trait_change('+input,b.input_change')
    def _set_input_change(self):
        print 'C: input change received'
        self.input_change = True





os.remove( 'c.pickle' )
print 'remove picke'
c = C()        
print 'initial c'
print c.b.a.value
#print 'changing value'
#c.b.a.value += 1
#print 'changed value'
#print c.b.a.value

print 'saving c'
import pickle
file_name = 'c.pickle'
file = open( file_name,'w')
pickle.dump( c.b, file )
file.close()

print 'changing value'
c.b.a.value += 1
print 'changed value'
print c.b.a.value

