
from scipy.linalg import *
from numpy import *

# ***********transpose a 4-dimension array:***************
# transpose(B_ijkl) ==> B_lkji : see example 1   
# 
# Transpose the selected columns is also possible 
#  
# Swap can work the same as transpose


# Example 1: 4-Dimension array

a = array(
      [
       [
        [    [ 0,  1],[5,6] ],
        [    [ 2,  3],[7,8] ]
       ],
        
       [
        [    [10, 11],[14,15] ],
        [    [12, 13],[16,17] ]
       ]
      ])

# dimension
print a.shape

print "a =", a

# swap columns 0 and 3 - s
c = a.swapaxes(0,3)
# swap columns 1 and 2
d = c.swapaxes(1,2)

print "d=", d

# transpose Matrix a 
b = a.transpose()
##default:
##b = a.transpose(3,2,1,0)

print "     "
print "     "
print "b =", b


#print "a[0,0,0,1]", a[0,0,0,1]
#print "b[0,1,0,0]", b[0,1,0,0]
#print "b[1,0,0,0]", b[1,0,0,0]


# Example 2: 3-Dimension array

x = arange(24).reshape(2,3,4)

print "x =", x

y = a.transpose()

print "y =", y

x = arange(24).reshape(4,3,2)




       
