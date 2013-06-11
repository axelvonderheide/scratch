from numpy import  *

def fn(x):
    return sin(x)

vsin = frompyfunc(fn, 1, 1)

a = 1
A = array([1, 0.3, 20.])
A2 = array([2., 2., 2.])
A_ = array([22., 33., 43.])


b = sin(a)
print 'sin(a) =\n' ,b
print '\n'

b = vsin(a)
print 'vsin(a) =\n' ,b
print '\n'

B = vsin(A)
print 'vsin(A) =\n' ,B
print '\n'


def myfunc(x, r):
   return r+(x**2)

print 'myfunc(20.,2) =\n',myfunc(20,2.)                                               # works fine
print '\n'

#myfunc(array([-2,2]))                                    # doesn't work, try it...



def vfunc(x, r):
#    vmyf = vectorize(myfunc, otypes=[float])              # declare the return type as float
    vmyf = frompyfunc(myfunc,2,1)
    return vmyf(x,r)

#vmyfunc = vectorize(myfunc, otypes=[float])              # declare the return type as float
#print 'vmyfunc(A, 2 ) =\n', vmyfunc(A,2) 
#print '\n'
#vfuncf = frompyfunc(myfunc,2,1)
#print 'vmyfuncf(A, 2 ) =\n', vmyfuncf(A,2) 
#print '\n'

# second argument of the vectorized function can be either a scalar or an array of the size of the first array
# argument and values equal to the scalar 
# in the first case the scalar argument is bradcasted to an array of the correct size 
print 'vfunc(A,2) = \n', vfunc(A,2)
print '\n'

print 'vfunc(A,A2) = \n', vfunc(A,A2)
print '\n'

# the second argument can be of course also an array with different values
# (vectorized function with two parameters)
print 'vfunc(A,A_) = \n', vfunc(A,A_)
print '\n'







## def myfunc(i,j):
##     return (i+1)*(j+4-i)

## a = fromfunction(myfunc, (4,4))

## a = arange(0,30)
## a.shape = (5,6)

## for i in xrange(a.shape[0]):
##     for j in xrange(a.shape[1]):
##         a[i,j]= (i+1)*(j+1)+(j+2)
##         print 'a[%g,%d]=%g'   % (i,j,a[i,j]),
##     print
