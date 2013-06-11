

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
            
    domain.coord_max = (3,0,0)
    domain.shape = (3,)
    ts.bcond_list =  [BCDof(var='u', dof = 0, value = 0.),
                      BCDof(var='u', dof = 1, link_dofs = [2], link_coeffs = [0.5] ),
                      BCDof(var='u', dof = 2, link_dofs = [3], link_coeffs = [1.] ),
                      BCDof(var='f', dof = 3, value = 1 ) ]
    # system solver
    u = tloop.eval()
    # expected solution
    print u
    u_ex = array([-0. ,  0.1 , 0.2 , 0.2],
                  dtype = float )
    difference = sqrt( norm( u-u_ex ) )
    print difference
    #self.assertAlmostEqual( difference, 0 )         
    #
    # '---------------------------------------------------------------'
    # 'Clamped bar with recursive constraints (displ at right end)'
    # 'u[1] = 0.5 * u[2], u[2] = 1.0 * u[3], u[3] = 1'
    ts.bcond_list =  [BCDof(var='u', dof = 0, value = 0.),
                      BCDof(var='u', dof = 1, link_dofs = [2], link_coeffs = [0.5] ),
                      BCDof(var='u', dof = 2, link_dofs = [3], link_coeffs = [1.] ),
                      BCDof(var='u', dof = 3, value = 1 ) ]
    u = tloop.eval()
    # expected solution
    print u
    u_ex = array([0. ,  0.5 , 1 ,  1], dtype = float )
    difference = sqrt( norm( u-u_ex ) )

    print difference
