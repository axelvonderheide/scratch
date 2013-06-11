
#from sys_matrix import SysSparseMtx, SysDenseMtx
from numpy import array, zeros, arange, array_equal, hstack, dot, sqrt
from scipy.linalg import norm

from mathkit.matrix_la.sys_mtx_assembly import SysMtxAssembly
from mathkit.matrix_la.coo_mtx import COOSparseMtx
from mathkit.matrix_la.dense_mtx import DenseMtx
import unittest

from ibvpy.api import \
    TStepper as TS, RTraceGraph, RTraceDomainListField, TLoop, \
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
from ibvpy.fets.fets2D.fets2D9q import FETS2D9Q
from ibvpy.fets.fets2D.fets2D4q8u import FETS2D4Q8U
from ibvpy.fets.fets1D.fets1D2l import FETS1D2L
from ibvpy.dots.dots_eval import DOTSEval
from ibvpy.dots.dots_list_eval import DOTSListEval

def ghost_bar( ):

    fets_eval = FETS1D2L( mats_eval = MATS1DElastic() ) 

    discr = ( 3, )
    # Discretization
    fe_domain1 = FEGrid( coord_min = (0,0,0),
                               coord_max = (2,0,0), 
                               shape   = discr,
                               inactive_elems = [0],
                               fets_eval = fets_eval )

    ts = TS( sdomain = fe_domain1,
             dof_resultants = True,
             bcond_list =  [ BCDof( var='u', value = 0., dof = 0 ),
                             BCDof( var='u', value = 0., dof = 3 ),
                             BCDof( var='f', value = 1., dof = 1 ) ],
            )

    # Add the time-loop control
    global tloop
    tloop = TLoop( tstepper = ts, tolerance = 1e-4, KMAX = 50,
                        tline  = TLine( min = 0.0,  step = 1.0, max = 1.0 ))

    print tloop.eval()
    
def L_shape( ):
    '''L-shaped domain constructed by deleting elements from the quadrangle'''
    
    mp = MATS2DElastic( E = 34.e3,
                        nu = 0.2 ) 

    fets_eval = FETS2D4Q(mats_eval = mp ) 

    discr = ( 3, 2 )
    # Discretization
    fe_domain1 = FEGrid( coord_min = (0,0,0),
                              coord_max = (2.,2.,0.), 
                              shape   = discr,
                              inactive_elems = [3,5],
                              fets_eval = fets_eval )

    ts = TS( sdomain = fe_domain1,
             dof_resultants = True,
             bcond_list =  [ BCDofGroup( var='u', value = 0., dims = [0,1],
                                         get_dof_method = fe_domain1.get_top_dofs ),
                             BCDofGroup( var='u', value = 0., dims = [0,1],
                                         get_dof_method = fe_domain1.get_left_dofs ),                                         
#                             BCDof( var='u', value = 0., dof = 20 ),
#                             BCDof( var='u', value = 0., dof = 21 ), 
#                             BCDof( var='u', value = 0., dof = 16 ), 
#                             BCDof( var='u', value = 0., dof = 17 ), 
                             BCDofGroup( var='f', value = -1, dims = [1],
                                         get_dof_method = fe_domain1.get_right_dofs ) 
                        ],
             rtrace_list =  [ RTraceDomainListField( name = 'Displacement', 
                                                     var = 'u', 
                                                     idx = 1, warp = False ),
                              RTraceDomainField(name = 'Stress' ,
                              var = 'sig_app', idx = 0,
                              record_on = 'update',
                              warp = True),
                              ]             
            )

    # Add the time-loop control
    global tloop
    tloop = TLoop( tstepper = ts, tolerance = 1e-4, KMAX = 50,
                        tline  = TLine( min = 0.0,  step = 1.0, max = 1.0 ))

    print tloop.eval()

    from ibvpy.plugins.ibvpy_app import IBVPyApp
    app = IBVPyApp( ibv_resource = tloop )
    app.main()      

def screwed_chess_board( ):
    '''Clamped bar 3 domains, each with 2 elems (displ at right end)
    [0]-[1]-[2] [3]-[4]-[5] [6]-[7]-[8]
    u[0] = 0, u[2] = u[3], u[5] = u[6], u[8] = 1'''
    
    mp = MATS2DElastic( E = 34.e3,
                        nu = 0.2 ) 

    fets_eval = FETS2D4Q(mats_eval = mp ) 
    #fets_eval = FETS2D4Q8U(mats_eval = mp ) 

    nx = 8
    ny = 8
    discr = ( nx, ny )
    inactive_elems = []
    for j in range( ny / 2 ):
        inactive_elems += [ i*2*ny+(j*2) for i in range( nx / 2 ) ] + \
                          [ (ny+1)+i*2*ny+(j*2) for i in range( nx / 2 ) ]
    load_dof = (ny+1) * 2 * (nx / 2) + ny
    # Discretization
    fe_domain1 = FEGrid( coord_min = (0,0,0),
                          coord_max = (2.,2.,0.), 
                          shape   = discr,
                          inactive_elems = inactive_elems,
                          fets_eval = fets_eval )

    ts = TS( sdomain = fe_domain1,
             dof_resultants = True,
             bcond_list =  [ BCDofGroup( var='u', value = 0., dims = [0,1],
                                         get_dof_method = fe_domain1.get_bottom_dofs ),
                             BCDofGroup( var='u', value = 0, dims = [0,1],
                                         get_dof_method = fe_domain1.get_top_dofs ), 
                             BCDofGroup( var='u', value = 0, dims = [0,1],
                                         get_dof_method = fe_domain1.get_left_dofs ), 
                             BCDofGroup( var='u', value = 0, dims = [0,1],
                                         get_dof_method = fe_domain1.get_right_dofs ), 
                             BCDof( var='f', value = 100., dof = load_dof ),
                             BCDof( var='f', value = 100., dof = load_dof+1 ),
                        ],
             rtrace_list =  [ RTraceDomainListField( name = 'Displacement', 
                                                     var = 'u', 
                                                     idx = 1, warp = False ),
                              RTraceDomainListField(name = 'Stress' ,
                              var = 'sig_app', idx = 0,
                              record_on = 'update',
                              warp = True),
                              ]             
            )

    # Add the time-loop control
    global tloop
    tloop = TLoop( tstepper = ts, tolerance = 1e-4, KMAX = 50,
                        tline  = TLine( min = 0.0,  step = 1.0, max = 1.0 ))

    tloop.eval()

    from ibvpy.plugins.ibvpy_app import IBVPyApp
    app = IBVPyApp( ibv_resource = tloop )
    app.main()      


if __name__ == '__main__':

#    ghost_bar()
#    L_shape()
    screwed_chess_board()
