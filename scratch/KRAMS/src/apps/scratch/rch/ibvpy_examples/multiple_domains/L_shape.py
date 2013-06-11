
#from sys_matrix import SysSparseMtx, SysDenseMtx
from numpy import array, zeros, arange, array_equal, hstack, dot, sqrt
from scipy.linalg import norm

from mathkit.matrix_la.sys_mtx_assembly import SysMtxAssembly
from mathkit.matrix_la.coo_mtx import COOSparseMtx
from mathkit.matrix_la.dense_mtx import DenseMtx
import unittest

from ibvpy.api import \
    TStepper as TS, RTraceGraph, RTraceDomainField, TLoop, \
    TLine, BCDof, BCDofGroup, IBVPSolve as IS

from ibvpy.mesh.fe_grid import FEGrid
from ibvpy.rtrace.rt_domain_list_field import \
     RTraceDomainListField 

from ibvpy.mats.mats1D.mats1D_elastic.mats1D_elastic import MATS1DElastic
from ibvpy.fets.fets1D.fets1D2l import FETS1D2L

from ibvpy.mats.mats2D.mats2D_elastic.mats2D_elastic import MATS2DElastic
from ibvpy.mats.mats2D.mats2D_cmdm.mats2D_cmdm import MA2DCompositeMicroplaneDamage
from ibvpy.mats.mats2D.mats2D_sdamage.mats2D_sdamage import MATS2DScalarDamage
from ibvpy.fets.fets2D.fets2D4q import FETS2D4Q
from ibvpy.mesh.fe_grid import FEGrid
from ibvpy.dots.dots_eval import DOTSEval
from ibvpy.dots.dots_list_eval import DOTSListEval


def L_shaped( ):
    '''Clamped bar 3 domains, each with 2 elems (displ at right end)
    [0]-[1]-[2] [3]-[4]-[5] [6]-[7]-[8]
    u[0] = 0, u[2] = u[3], u[5] = u[6], u[8] = 1'''
    
    mp = MATS2DScalarDamage(E = 34.e3,
                                   nu = 0.2,
                                   epsilon_0 = 59.e-6,
                                   epsilon_f = 3.2e-3,
                                   #epsilon_f = 3.2e-1,
                                   #stiffness  = "algorithmic",
                                   strain_norm_type = 'Mises')
                                                   
    mp = MA2DCompositeMicroplaneDamage()
#    mp = MATS2DElastic( E = 34.e3,
#                        nu = 0.2 ) 

    fets_eval = FETS2D4Q(mats_eval = mp ) 
    
    discr = ( 10, 10 )
    # Discretization
    fe_grid1 = FEGrid( coord_min = (0,0,0),
                               coord_max = (1.,1.,0.), 
                               shape   = discr,
                               fets_eval = fets_eval )

    fe_grid2 = FEGrid( coord_min = (0.,1.,0),
                               coord_max = (1.,2.,0.), 
                               shape   = discr,
                               fets_eval = fets_eval )
    
    fe_grid3 = FEGrid( coord_min = (1.,1.,0),
                               coord_max = (2.,2.,0.), 
                               shape   = discr,
                               fets_eval = fets_eval )

    ts = TS( sdomain = [ fe_grid1, fe_grid2, fe_grid3 ],
             dof_resultants = True,
             bcond_list =  [ BCDofGroup( var='u', value = 0., dims = [0,1],
                                         get_dof_method = fe_grid1.get_bottom_dofs ),
                             BCDofGroup( var='u', value = 0., dims = [0,1],
                                         get_dof_method = fe_grid3.get_left_dofs,
                                         get_link_dof_method = fe_grid2.get_right_dofs,
                                         link_coeffs = [1.] ),                                                                        
                             BCDofGroup( var='u', value = 0., dims = [0,1],
                                         get_dof_method = fe_grid2.get_bottom_dofs,
                                         get_link_dof_method = fe_grid1.get_top_dofs,
                                         link_coeffs = [1.] ),
                             BCDofGroup( var='u', value = 0.0004, dims = [1],
                                         get_dof_method = fe_grid3.get_right_dofs ) ],
             rtrace_list =  [ RTraceDomainListField( name = 'Displacement', 
                                                     var = 'u', 
                                                     idx = 1 ),
#                              RTraceDomainListField(name = 'Damage' ,
#                              var = 'omega', idx = 0,
#                              record_on = 'update',
#                              warp = False ),
                              RTraceDomainListField(name = 'Stress' ,
                              var = 'sig_app', idx = 0,
                              record_on = 'update',
                              warp = True),
                              RTraceDomainListField(name = 'Phi in principle damage coords' ,
                              var = 'phi_pdc', idx = 0,
                              record_on = 'update',
                              warp = True),
                              RTraceDomainListField(name = 'Fracture energy' ,
                              var = 'fracture_energy', idx = 0,
                              record_on = 'update',
                              warp = True),
#                              RTraceDomainListField(name = 'Strain' ,
#                              var = 'eps_app', idx = 0,
#                              record_on = 'update',
#                              warp = False), 
                              ]             
            )

    # Add the time-loop control
    global tloop
    tloop = TLoop( tstepper = ts, tolerance = 1e-2, KMAX = 200,
                        tline  = TLine( min = 0.0,  step = 1.0, max = 1.0 ))

    tloop.eval()
#
    from ibvpy.plugins.ibvpy_app import IBVPyApp
    app = IBVPyApp( ibv_resource = tloop )
    app.main()      


if __name__ == '__main__':

    L_shaped()