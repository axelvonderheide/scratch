
from scipy.linalg import *
from numpy import *

#x = array([[1,2,3],[4,5,6]])
#y = array([[1,2],[3,4],[5,6]])
# 
#print dot(x,y)
#
#
#x = array([1,2,3])
#y = array([10,20,30])
#print inner(x,y)    
#
x = array([1,2,3])
y = array([10,20,30])
print outer(x,y) 

#
#a = arange(60.).reshape(3,4,5)
#b = arange(24.).reshape(4,3,2)
#c = tensordot(a,b, axes=([1,0],[0,1]))
#
#print c
#
#
#print arange(3)
#
#
#m = array([2,2,2,2])
#n = array([1,2,3])



#a = arange(4.).reshape(2,2)
#b = arange(4.).reshape(2,2)

#a = array([[1,1],[2,3]])
#b = array([[1,1],[4,5]])
#
#c = tensordot(a,b, axes=([1,0],[1,0]))
#
#print c



#a = array(
#      [
#       [
#        
#        [    [ 0,  1],[5,6] ],
#        [    [ 2,  3],[7,8] ]
#        
#        
#       ],
#        
#       [
#        
#        [    [10, 11],[14,15] ],
#        [    [12, 13],[16,17] ]
#        
#       ]
#       
#      ])
#
#print a.shape
#print "a =", a
#
#c = a.swapaxes(0,3)
#
#d = c.swapaxes(1,2)
#
#print "d=", d
#
#
#b = a.transpose()
###default:
##b = a.transpose(3,2,1,0)
#
#print "     "
#print "     "
#print "b =", b
#
#
#print "a[0,0,0,1]", a[0,0,0,1]
#print "b[1,0,0,0]", b[1,0,0,0]
#print "b[0,1,0,0]", b[0,1,0,0]




#a = arange(24).reshape(2,3,4)
#
#print "a =", a
#
#b = a.transpose()
#
#print "b =",b
#
#a = arange(24).reshape(4,3,2)



