'''
Created on Mar 31, 2009

@author: jakub
'''

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

from ibvpy.mesh.fe_grid_domain import FEGridDomain
from ibvpy.rtrace.rt_domain_list_field import \
     RTraceDomainListField 

from ibvpy.mats.mats1D.mats1D_elastic.mats1D_elastic import MATS1DElastic
from ibvpy.fets.fets1D.fets1D2l import FETS1D2L

from ibvpy.mats.mats2D.mats2D_elastic.mats2D_elastic import MATS2DElastic
from ibvpy.mats.mats2D.mats2D_cmdm.mats2D_cmdm import MA2DCompositeMicroplaneDamage
from ibvpy.mats.mats2D.mats2D_sdamage.mats2D_sdamage import MATS2DScalarDamage
from ibvpy.fets.fets2D.fets2D4q import FETS2D4Q
from ibvpy.mesh.fe_domain_list import FEDomainList
from ibvpy.dots.dots_eval import DOTSEval
from ibvpy.dots.dots_list_eval import DOTSListEval

def xtest_notched( ):
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
                                                   
#    mp = MATS2DElastic( E = 34.e3,
#                        nu = 0.2 ) 
    fets_eval = FETS2D4Q(mats_eval = mp ) 

    discr = ( 1, 1 )
    # Discretization
    fe_domain1 = FEGridDomain( coord_min = (0,0,0),
                               coord_max = (0.203,0.082,0), 
                               shape   = discr,
                               fets_eval = fets_eval )

    fe_domain2 = FEGridDomain( coord_min = (0,0.082,0),
                               coord_max = (0.203,0.306,0), 
                               shape   = discr,
                               fets_eval = fets_eval )
    
    fe_domain3 = FEGridDomain( coord_min = (0.203,0,0),
                               coord_max = (0.6,0.082,0), 
                               shape   = (20,2),
                               fets_eval = fets_eval )
    
    fe_domain4 = FEGridDomain( coord_min = (0.203,0.082,0),
                               coord_max = (0.6,0.306,0), 
                               shape   = (20,2),
                               fets_eval = fets_eval )
    
    fe_domain5 = FEGridDomain( coord_min = (0.6,0,0),
                               coord_max = (0.661,0.082,0), 
                               shape   = (20,2),
                               fets_eval = fets_eval )

    fe_domain6 = FEGridDomain( coord_min = (0.6,0.082,0),
                               coord_max = (0.661,0.306,0), 
                               shape   = (20,2),
                               fets_eval = fets_eval )
    
    fe_domain7 = FEGridDomain( coord_min = (0.661,0,0),
                               coord_max = (0.722,0.082,0), 
                               shape   = (20,2),
                               fets_eval = fets_eval )
    
    fe_domain8 = FEGridDomain( coord_min = (0.661,0.082,0),
                               coord_max = (0.722,0.306,0), 
                               shape   = (20,2),
                               fets_eval = fets_eval )
    
    fe_domain9 = FEGridDomain( coord_min = (0.722,0,0),
                               coord_max = (1.119,0.082,0), 
                               shape   = (20,2),
                               fets_eval = fets_eval )
    
    fe_domain10 = FEGridDomain( coord_min = (0.722,0.082,0),
                               coord_max = (1.119,0.306,0), 
                               shape   = (20,2),
                               fets_eval = fets_eval )
    
    fe_domain11 = FEGridDomain( coord_min = (1.119,0,0),
                               coord_max = (1.322,0.082,0), 
                               shape   = (20,2),
                               fets_eval = fets_eval )
    
    fe_domain12 = FEGridDomain( coord_min = (1.119,0.082,0),
                               coord_max = (1.322,0.306,0), 
                               shape   = (20,2),
                               fets_eval = fets_eval )

    fe_domain  = FEDomainList( subdomains = [ fe_domain1, fe_domain2, fe_domain3, fe_domain4,\
                                              fe_domain5, fe_domain6, fe_domain7, fe_domain8,\
                                              fe_domain9, fe_domain10, fe_domain11, fe_domain12])

    ts = TS( sdomain = fe_domain,
             dof_resultants = True,
             bcond_list =  [ BCDofGroup( var='u', value = 0., dims = [0,1],
                                         get_dof_method = fe_domain3.get_bottom_right_dofs ),
                             BCDofGroup( var='u', value = 0., dims = [1],
                                         get_dof_method = fe_domain9.get_bottom_right_dofs ),
                                         
                             BCDofGroup( var='u', value = 0., dims = [0,1],
                                         get_dof_method = fe_domain2.get_bottom_dofs,
                                         get_link_dof_method = fe_domain1.get_top_dofs,
                                         link_coeffs = [1.] ),                                                                        
                             BCDofGroup( var='u', value = 0., dims = [0,1],
                                         get_link_dof_method = fe_domain4.get_bottom_dofs,
                                         get_dof_method = fe_domain3.get_top_dofs,
                                         link_coeffs = [1.] ),   
                            BCDofGroup( var='u', value = 0., dims = [0,1],
                                         get_dof_method = fe_domain6.get_bottom_dofs,
                                         get_link_dof_method = fe_domain5.get_top_dofs,
                                         link_coeffs = [1.] ),                                                                        
                             BCDofGroup( var='u', value = 0., dims = [0,1],
                                         get_dof_method = fe_domain8.get_bottom_dofs,
                                         get_link_dof_method = fe_domain7.get_top_dofs,
                                         link_coeffs = [1.] ),
                            BCDofGroup( var='u', value = 0., dims = [0,1],
                                         get_dof_method = fe_domain10.get_bottom_dofs,
                                         get_link_dof_method = fe_domain9.get_top_dofs,
                                         link_coeffs = [1.] ),                                                                        
                             BCDofGroup( var='u', value = 0., dims = [0,1],
                                         get_dof_method = fe_domain12.get_bottom_dofs,
                                         get_link_dof_method = fe_domain11.get_top_dofs,
                                         link_coeffs = [1.] ),                           
                                         
                                         
                            BCDofGroup( var='u', value = 0., dims = [0,1],
                                         get_link_dof_method = fe_domain3.get_left_dofs,
                                         get_dof_method = fe_domain1.get_right_dofs,
                                         link_coeffs = [1.] ),                                                                        
                             BCDofGroup( var='u', value = 0., dims = [0,1],
                                         get_dof_method = fe_domain4.get_left_dofs,
                                         get_link_dof_method = fe_domain2.get_right_dofs,
                                         link_coeffs = [1.] ),
                                         
                            BCDofGroup( var='u', value = 0., dims = [0,1],
                                         get_dof_method = fe_domain5.get_left_dofs,
                                         get_link_dof_method = fe_domain3.get_right_dofs,
                                         link_coeffs = [1.] ),                                                                        
                             BCDofGroup( var='u', value = 0., dims = [0,1],
                                         get_link_dof_method = fe_domain6.get_left_dofs,
                                         get_dof_method = fe_domain4.get_right_dofs,
                                         link_coeffs = [1.] ),
#                                        
#                            BCDofGroup( var='u', value = 0., dims = [0,1],
#                                         get_dof_method = fe_domain5.get_left_dofs,
#                                         get_link_dof_method = fe_domain3.get_right_dofs,
#                                         link_coeffs = [1.] ),                                                                        
#                             BCDofGroup( var='u', value = 0., dims = [0,1],
#                                         get_dof_method = fe_domain6.get_left_dofs,
#                                         get_link_dof_method = fe_domain4.get_right_dofs,
#                                         link_coeffs = [1.] ),
#                            
#                            BCDofGroup( var='u', value = 0., dims = [0,1],
#                                         get_dof_method = fe_domain8.get_left_dofs,
#                                         get_link_dof_method = fe_domain6.get_right_dofs,
#                                         link_coeffs = [1.] ),   
#                                         
#                            BCDofGroup( var='u', value = 0., dims = [0,1],
#                                         get_dof_method = fe_domain9.get_left_dofs,
#                                         get_link_dof_method = fe_domain7.get_right_dofs,
#                                         link_coeffs = [1.] ),                                                                        
#                             BCDofGroup( var='u', value = 0., dims = [0,1],
#                                         get_dof_method = fe_domain10.get_left_dofs,
#                                         get_link_dof_method = fe_domain8.get_right_dofs,
#                                         link_coeffs = [1.] ),
#                                         
#                            BCDofGroup( var='u', value = 0., dims = [0,1],
#                                         get_dof_method = fe_domain11.get_left_dofs,
#                                         get_link_dof_method = fe_domain9.get_right_dofs,
#                                         link_coeffs = [1.] ),                                                                        
#                             BCDofGroup( var='u', value = 0., dims = [0,1],
#                                         get_dof_method = fe_domain12.get_left_dofs,
#                                         get_link_dof_method = fe_domain10.get_right_dofs,
#                                         link_coeffs = [1.] ),             
#                                         
#                            BCDofGroup( var='f', value = -1., dims = [1],
#                                         get_dof_method = fe_domain8.get_top_right_dofs)
                                         
                              ],
             rtrace_list =  [ 
#                             RTraceDomainListField( name = 'Displacement', 
#                                                     var = 'u', 
#                                                     idx = 1 ),
                              RTraceDomainListField(name = 'Damage' ,
                              var = 'omega', idx = 0,
                              record_on = 'update',
                              warp = True),
                              RTraceDomainListField(name = 'Stress' ,
                              var = 'sig_app', idx = 0,
                              record_on = 'update'),
#                              RTraceDomainListField(name = 'Strain' ,
#                              var = 'eps_app', idx = 0,
#                              record_on = 'update',
#                              warp = False), 
                              ]             
            )

    # Add the time-loop control
    global tloop
    tloop = TLoop( tstepper = ts, tolerance = 1e-4, KMAX = 50,
                        tline  = TLine( min = 0.0,  step = 1., max = 1.0 ))

    tloop.eval()
#    import cProfile
#    cProfile.run('tloop.eval()', 'tloop_prof' )
#    
#    import pstats
#    p = pstats.Stats('tloop_prof')
#    p.strip_dirs()
#    print 'cumulative'
#    p.sort_stats('cumulative').print_stats(20)
#    print 'time'
#    p.sort_stats('time').print_stats(20)    

    from ibvpy.plugins.ibvpy_app import IBVPyApp
    app = IBVPyApp( ibv_resource = tloop )
    app.main()      


if __name__ == '__main__':

    xtest_notched()