'''
exTreated on exTMpr 30, 2010

@author: alexander
'''

from enthought.traits.api import \
    HasTraits, Instance, Int, on_trait_change, \
    Event, DelegatesTo

import pickle

import os

from os.path import \
    exists

#--------------------------------------
# Without the delegated attribute 'rho' 
# of class B the input change is propagated 
# and received correctly as expected.
#
# --- 1. question: ---
# why doesn't this work if the delegated 
# attribute 'rho' is defined in 'B'? 
# 'rho' itself is not even modified ?
# 
# --- 2. question: ---
# why does the input change is received if
# the delegated attribute is named 'rrho'
# instead of 'rho' (or '_value' instead of 'value')?
# How does the name effect the behavior of 
# 'input_change'in 'B'?
#
# --- 3. question: ---
# why does the input change is received 
# correctly even with a delegated attribute
# named 'rho' if the class 'B' is read in 
# from a pickled file?
# 
# --- 4. question: ---
# If 'listenable = False' is specified in DelegatesTo
# the changes are propagated correctly!
#
#--------------------------------------

# CCS
class A( HasTraits ):

    # not changed/unused but delegated to B
    rho = Int(3, input = True)
    rrho = Int(3, input = True)
    
    value  = Int( 1, input = True )

    input_change = Event
    @on_trait_change('+input')
    def _set_input_change(self):
        print 'A: input change propagated !!!'
        self.input_change = True

# CTT
class B( HasTraits ):
    
    a = Instance( A, input = True )
    def _a_default(self):
        return A()

    input_change = Event
    @on_trait_change('a.input_change')
    def _set_input_change(self):
        print 'B: input change received !!!'
        self.input_change = True

    # NOTE: if 'rho' is defined no input change is received! Why?
    rho = DelegatesTo('a')

    # attribute named 'rrho' has no effect
    rrho = DelegatesTo('a')
    
    # with option 'listenable = False' notification of the 
    # delegeted attribute is explicity turned off and 
    # 'on_trait_change' works correctly even with name 'rho'
#    rho = DelegatesTo('a', listenable=False )



#---------------------------------------------------
# delete pickle file if it exists:
#---------------------------------------------------
print '\n' 
print 'XXX delete pickle file if it exists '

if os.path.exists( 'b.pickle' ):
    print 'pickle file removed'
    os.remove( 'b.pickle' )
else:
    print 'pickle file does not exist'
    
#---------------------------------------------------
# Construct b and show default attribute:
#---------------------------------------------------
print '\n' 
print '--- STATE 1: ---' 
print 'XXX Construct b and show default value'

b = B()        

print 'b.a.value = ', b.a.value

#---------------------------------------------------
# Save b as pickle:
#---------------------------------------------------
print 'XXX Save b as pickle'

print 'saving b'
file_name = 'b.pickle'
file = open( file_name,'w')
pickle.dump( b, file )
file.close()

#---------------------------------------------------
# change attribute:
#---------------------------------------------------
print 'XXX change value of b.a.value'

b.a.value += 1
print 'b.a.value = ', b.a.value


#---------------------------------------------------
# Save b as pickle:
#---------------------------------------------------
print 'XXX Save b as pickle'

print 'saving b'
file_name = 'b.pickle'
file = open( file_name,'w')
pickle.dump( b, file )
file.close()

#---------------------------------------------------
# load b from pickle:
#---------------------------------------------------
print '\n' 
print '--- STATE 2: ---' 
print 'XXX load b from pickle and show default value'

file = open( 'b.pickle', 'r' )
b = pickle.load( file )

print 'b.a.value = ', b.a.value

#---------------------------------------------------
# change attribute:
#---------------------------------------------------
print 'XXX change value of b.a.value'

b.a.value += 1
print 'b.a.value = ', b.a.value

#---------------------------------------------------
# clean up
#---------------------------------------------------
#if os.path.exists( 'b.pickle' ):
#    os.remove( 'b.pickle' )
