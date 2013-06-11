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
# Created on Aug 18, 2009 by: rch

from ibvpy.fets.fets_eval import FETSEval

from numpy import \
     zeros, dot, hstack, identity

from scipy.linalg import \
     inv     

#-----------------------------------------------------------------------------------
# FETS3D - Base class for 3D elements 
#-----------------------------------------------------------------------------------

n_nodes = 2
n_dim = 2

if n_dim == 2:
    n_strain = 3
if n_dim == 3:
    n_strain = 6

dNx_mtx = arange(n_dim*n_nodes).reshape(n_dim,n_nodes)
print 'dNx_mtx', dNx_mtx
    
Bx_mtx = zeros( (n_strain, n_nodes * n_dim ), dtype = 'float_' )      
Bx_mtx[[0,4,5], ::3]  = dNx_mtx[[0,2,1]]
Bx_mtx[[1,3,5], 1::3] = dNx_mtx[[1,2,0]]
Bx_mtx[[2,3,4], 2::3] = dNx_mtx[[2,1,0]]
return Bx_mtx


#class FETS3D(FETSEval):
#    '''Base class for 3D elements.
#    '''
#    # Dimensional mapping
#    dim_slice = slice(0, 3)
#    
#    def get_B_mtx( self, r_pnt, X_mtx ):
#        '''
#        Return the kinematic matrix
#        @param r_pnt:local coordinates
#        @param X_mtx:global coordinates
#        '''
#        J_mtx = self.get_J_mtx(r_pnt,X_mtx)
#        dNr_mtx = self.get_dNr_mtx( r_pnt )
#        dNx_mtx = dot( inv( J_mtx ), dNr_mtx  )
#        n_nodes = len( self.dof_r )
#        Bx_mtx = zeros( (6, n_nodes * 3 ), dtype = 'float_' )      
#        Bx_mtx[[0,4,5], ::3]  = dNx_mtx[[0,2,1]]
#        Bx_mtx[[1,3,5], 1::3] = dNx_mtx[[1,2,0]]
#        Bx_mtx[[2,3,4], 2::3] = dNx_mtx[[2,1,0]]
#        return Bx_mtx
    