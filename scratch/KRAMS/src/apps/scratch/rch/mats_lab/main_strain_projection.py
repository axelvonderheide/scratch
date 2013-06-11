#-------------------------------------------------------------------------------
#
# Copyright (c) 2009, IMB, RWTH Aachen.
# All rights reserved.
#
# This software is provided without warranty under the terms of the BSD
# license included in simvisage/LICENSE.txt and may be redistributed only
# under the conditions described in the aforementioned license.  The license
# is also available online at http://www.simvisage.com/licenses/BSD.txt
#
# Thanks for using Simvisage open source!
#
# Created on Sep 3, 2009 by: rch

from ibvpy.mats.mats3D.mats3D_tensor import map3d_eps_eng_to_mtx
from numpy import min, max, dot, array, argmax
from scipy.linalg import eigh
from math import sqrt

eps_app_eng = array( [0.1, 0.2, 0.8
                      , 1.0, 0.0, 0.0],dtype = 'float_' )
X_mtx = array( [[-1, -1, -1],
                [-1, -1,  1],
                [-1,  1, -1],
                [-1,  1,  1],
                [ 1, -1, -1],
                [ 1, -1,  1],
                [ 1,  1, -1],
                [ 1,  1,  1]], dtype = 'float_',
                 )

# first principle strain unit vector
#        
eps_mtx = map3d_eps_eng_to_mtx( eps_app_eng )
print 'eps_mtx',eps_mtx

eigval, eigvec = eigh( eps_mtx )
print 'eigh',eigval
eig1 = argmax(eigval)

eps_one = eigvec[eig1]
print 'eps_one', eps_one
proj = dot( X_mtx, eps_one )
h = max( proj )-min( proj )
print 'diagonal',2* sqrt(2)

print h