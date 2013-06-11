from scipy.linalg import *
import numpy
import scipy.linalg

from numpy import array, ones, zeros, outer, inner, transpose, dot, frompyfunc, \
     fabs, sqrt, linspace, vdot, identity, tensordot, sin as nsin, meshgrid, float_, ix_, \
     vstack, hstack, sqrt as arr_sqrt

from scipy.linalg import eig, inv


phi_mtx = array([[1,6,5],[6,2,4],[5,4,3]])

phi_eig_val, phi_eig_mtx = eig( phi_mtx )
phi_eig = array([ pe.real for pe in phi_eig_val] )

# use 'dot' 
phi_pdc_mtx = dot(dot(transpose(phi_eig_mtx),phi_mtx),phi_eig_mtx)

# use 'tensordot' (false):
#phi_pdc_mtx_ten = tensordot( tensordot( phi_mtx, phi_eig_mtx, [[0],[0]] ), phi_eig_mtx, [[1],[0]] )

# use 'tensordot (correct):
phi_pdc_mtx_ten = tensordot( tensordot( phi_eig_mtx, phi_mtx, [[0],[0]] ), phi_eig_mtx, [[1],[0]] )


## # verify the transformation
## #phi_pdc_mtx = tensordot( tensordot( phi_eig_mtx, phi_mtx, [[0],[0]] ), phi_eig_mtx, [[1],[0]] )
##     phi_pdc_mtx = tensordot( tensordot( phi_mtx, phi_eig_mtx, [[0],[0]] ), phi_eig_mtx, [[1],[0]] )
## w_mtx = arr_sqrt( phi_pdc_mtx )

n_d = 3
phi_pdc_mtx_v = identity(n_d)
for i in range(0,n_d):
    phi_pdc_mtx_v[i,i]=phi_eig[i]

w_mtx_pdc = arr_sqrt( phi_pdc_mtx_v )
    
w_mtx = dot(dot(phi_eig_mtx,w_mtx_pdc),transpose(phi_eig_mtx))



print "\n"

print "phi_mtx = ", phi_mtx ,"\n"
print "phi_eig_val = ", phi_eig_val ,"\n"
print "phi_eig = ", phi_eig ,"\n"
print "phi_eig_mtx = ", phi_eig_mtx ,"\n"
print "use dot: phi_pdc_mtx = ", phi_pdc_mtx, "\n"
print "use tensordot (correct order of arguments): phi_pdc_mtx_ten = ", phi_pdc_mtx_ten, "\n"
print "construct phi_pdc based on the eigen-values: phi_pdc_mtx_v = ", phi_pdc_mtx_v, "\n"
print "tensorial squareroot: w_mtx_pdc = ", w_mtx_pdc, "\n"
print "w_mtx = ", w_mtx, "\n"




phi_pdc_mtx = dot(dot(transpose(phi_eig_mtx),phi_mtx),phi_eig_mtx)
