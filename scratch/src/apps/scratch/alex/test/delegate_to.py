'''
Created on May 3, 2010

@author: alexander
'''
from enthought.traits.api \
    import DelegatesTo, HasTraits, Instance, Str

class Parent(HasTraits):
    first_name = Str
    last_name  = Str

class Child(HasTraits):
    first_name = Str
    last_name  = DelegatesTo('father', listenable = False)
    father     = Instance(Parent)
    mother     = Instance(Parent)

""
tony  = Parent(first_name='Anthony', 
               last_name='Jones')

alice = Parent(first_name='Alice', 
               last_name='Smith')

sally = Child( first_name='Sally', 
               father=tony, 
               mother=alice)

print sally.last_name
print tony.last_name

sally.last_name = 'Cooper' # Updates delegatee

print sally.last_name
print tony.last_name

tony.last_name = 'AAA' # Updates delegate

print sally.last_name
print tony.last_name

