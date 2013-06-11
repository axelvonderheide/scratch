

import sys
print sys.path
#from sys_matrix import SysSparseMtx, SysDenseMtx
from numpy import array, zeros, arange, array_equal, hstack, dot, sqrt
from scipy.linalg import norm

from mathkit.matrix_la.sys_mtx_assembly import SysMtxAssembly
from mathkit.matrix_la.coo_mtx import COOSparseMtx
from mathkit.matrix_la.dense_mtx import DenseMtx
import unittest

from ibvpy.api import \
    TStepper as TS, RTraceGraph, RTraceDomainField, TLoop, \
    TLine, BCDof, IBVPSolve as IS, DOTSEval
from ibvpy.mats.mats1D.mats1D_elastic.mats1D_elastic import MATS1DElastic

from ibvpy.bcond.bc_slice import BCSlice
from ibvpy.mesh.fe_grid import FEGrid
from ibvpy.fets.fets1D.fets1D2l import FETS1D2L

if __name__ == '__main__':
    
    fets_eval = FETS1D2L(mats_eval = MATS1DElastic(E=10., A=1.))        
    
    # Discretization
    domain = FEGrid( coord_max = (10.,0.,0.), 
                     shape   = (1,),
                     fets_eval = fets_eval )
    
    ts = TS( sdomain = domain,
                  dof_resultants = True
                        )
    tloop = TLoop( tstepper = ts,
                        tline  = TLine( min = 0.0,  step = 1, max = 1.0 ))
    
        
    domain.coord_max = (10,0,0)
    domain.shape = (10,)
    bc_left  = BCSlice( var = 'u', value = 0., slice = domain[ 0, 0] )
    bc_right = BCSlice( var = 'u', value = 1., slice = domain[1:,:] )
    ts.bcond_list =  [ bc_left, bc_right ]

#    ts.bcond_list =  [BCDof(var='u', dof = 0, value = 0.),
#                      BCDof(var='u', dof = 1, link_dofs = [2], link_coeffs = [0.5] ),
#                      BCDof(var='u', dof = 2, link_dofs = [3], link_coeffs = [1.] ),
#                      BCDof(var='u', dof = 3, value = 1. ) ]
    
    ts.rtrace_list = [ RTraceGraph(name = 'Fi,right over u_right (iteration)' ,
                                   var_y = 'F_int', idx_y = 10,
                                   var_x = 'U_k',   idx_x = 10) ]
    
    u = tloop.eval()
    print 'u',u
    # expected solution
    u_ex = array([ 0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1. ],
                 dtype = float )
    print 'u_ex',u_ex
    difference = sqrt( norm( u-u_ex ) )
    print 'difference', difference
    
    # compare the reaction at the left end
    F = ts.F_int[0]
    print 'F_int should be -1', F
