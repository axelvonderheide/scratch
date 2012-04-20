from numpy import allclose, arange, eye, linalg, ones, ix_, array
from scipy import linsolve, sparse

#Asp = sparse.lil_matrix((1000,1000))
Asp = sparse.dok_matrix((1000,1000))
Bsp = sparse.dok_matrix((10,10))
Asp.setdiag(ones(1000))
Afull = eye(1000)

print "start"
a = array([1,3,5])
b = array([12,15,20])

Afull[ix_(a),ix_(a)]=1.
Asp[a,b] = 1.

Bsp = Asp + Asp

#print Asp[ix_(a),ix_(b)]
Asp[2,:] = 3

Bsp = Asp + Asp
#Asp[ix_(a),5]=1.
print "after"
b = arange(1,1001)
Asp = Asp.astype('d')
xsp = linsolve.spsolve(Asp,b)
print xsp
print "finished"

Asp = sparse.dok_matrix((4,4))
Asp[2,:] = 5
Bsp = Asp + Asp

Asp[3,2] = 8

print 'col 2', Asp.getcol(2)
Bsp = Bsp + Asp

print 'Asp'
print Bsp[3,1]

print Bsp[:,2]