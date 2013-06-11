from numpy import *

A = arange(4)
print 'A', A

B = 1./A
print 'B', B

At = transpose(A)
Bt = transpose(B)
print 'At', At
print 'Bt', Bt


print '\n vstack(A,B)', transpose(vstack((A,B)))
