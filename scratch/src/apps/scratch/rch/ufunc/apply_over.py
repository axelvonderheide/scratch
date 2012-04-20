
# repeatedly apply the same function - 
from numpy import *
a = arange(24).reshape(2,3,4)         # a has 3 axes: 0,1 and 2
a

print 'a'
print a

def fn( a, b ):
    print 'a',a
    print 'b',b
    return a

app = apply_over_axes( fn, a, [1,2] )
print app

app = apply_over_axes( sum, a, [1,2] )

print 'sum',app

print sum( sum(a,1), 1)