
from numpy import *
from scipy.linalg import *

import numpy
import scipy.linalg


a = array([1,2])

print '\n '
print 'norm af a 1d array \n'
print 'XXX', norm(a)



A = numpy.array([[1,6,5],[6,2,4],[5,14,3]])
det_A = scipy.linalg.det(A)
inv_A = scipy.linalg.inv(A)

V,E = scipy.linalg.eig(A)
V_n,E_n = numpy.linalg.eig(A)

Id = numpy.identity(3)
A_1 = A-V[0]*Id
det_A_1 = det(A_1)

print "A = ", A, "\n"
print "det(A) = ", det_A, "\n"
print "inv(A) = ", inv_A, "\n"
print "Eigenvalues V = ", V ,"\n"
print "Eigenvectors E = ", E ,"\n"
print "use numpy insteadt of scipy: E = ", E_n ,"\n"
print "Identity matrix = ", Id ,"\n"
print "A_1 = ", A_1, "\n"
print "det(A_1)", det_A_1 
print "transpose(inv(A))", transpose(inv(A)) 
print "inv(transpose(A))", inv(transpose(A)) 

print '\n '
print 'test the numpy multiply-functions: \n'
a = array([1,2])
b = array([3,4])
A = array([[1,2],[0,3]])
B = array([[3,4],[2,3]])



#A_list = array([A,A,A])

#
#print 'A_list.sum', A_list.sum(0)

#print 'dot(A,transpose(A))', dot(A,transpose(A))
#
#print 'A*transpose(A)', A*transpose(A)

#print "a = ", a.shape
print "a = ", a
print "b = ", b, "\n"
#print "A = \n", A
#print "B = \n", B, "\n"

print "a*b = ", a*b


#print "inner(A,b) = ", inner(A,b)
#print "inner(a,b) = ", inner(a,b)
#print "inner(A,B) = ", inner(A,B), "\n"

#print "outer(a,b) = ", outer(a,b) 
#print "outer(A,B) = ", outer(A,B), "\n"
#
print "dot(A,b) = ", dot(A,b)
print "dot(a,b) = ", dot(a,b)
print "dot(a,bT) = ", dot(transpose(b),a)
#
#print "b = ", b
#print "bT = ", transpose(b)

print "dot(A,B) = ", dot(A,B), "\n"

print "A*B = ", A*B, "\n"



#print "vdot(a,b) = ", vdot(a,b)
#print "vdot(A,B) = ", vdot(A,B), "\n"
#
#
#print "inv(A) = ", inv(A), "\n"
#print "transpose(A) = ", transpose(A), "\n"
#
#
#
#
#
print "vdot(a,b) = ", vdot(a,b)
print "vdot(A,B) = ", vdot(A,B), "\n"

#print "tensordot(A,B,[[0,0],[1,1]]) = ", tensordot(a,b,[[0,1],[1,0]])
print "tensordot(A,B,[[0],[1]]) = ", tensordot(A,B,[[0],[1]]), "\n"
print "tensordot(A,B,[[1],[0]]) = ", tensordot(A,B,[[1],[0]]), "\n"

print "tensordot(a,B,[[0],[1]]) = ", tensordot(a,B,[[0],[1]]), "\n"
