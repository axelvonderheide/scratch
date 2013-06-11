#from sys_matrix import SysSparseMtx, SysDenseMtx
from numpy import array, zeros, arange, array_equal, hstack, dot, sqrt
from scipy.linalg import norm

from mathkit.matrix_la.sys_mtx_assembly import SysMtxAssembly
from mathkit.matrix_la.coo_mtx import COOSparseMtx
from mathkit.matrix_la.dense_mtx import DenseMtx
import unittest

from ibvpy.api import \
    TStepper as TS, RTraceGraph, RTraceDomainField, TLoop, \
    TLine, BCDof, BCDofGroup, IBVPSolve as IS, DOTSEval
from ibvpy.rtrace.rt_domain_list_field import RTraceDomainListField

from ibvpy.mesh.fe_grid import FEGrid
from ibvpy.mesh.fe_domain_list import FEDomainList
from ibvpy.dots.dots_list_eval import DOTSListEval

from ibvpy.mats.mats1D.mats1D_elastic.mats1D_elastic import MATS1DElastic
from ibvpy.fets.fets1D.fets1D2l import FETS1D2L

from ibvpy.mats.mats2D.mats2D_elastic.mats2D_elastic import MATS2DElastic
from ibvpy.mats.mats2D.mats2D_cmdm.mats2D_cmdm import MA2DCompositeMicroplaneDamage
from ibvpy.mats.mats2D.mats2D_sdamage.mats2D_sdamage import MATS2DScalarDamage
from ibvpy.fets.fets2D.fets2D4q import FETS2D4Q
from ibvpy.mats.mats1D5.mats1D5bond import MATS1D5Bond
from mathkit.mfn.mfn_line.mfn_line import MFnLineArray
from ibvpy.fets.fets1D5.fets1D52b4uLRH import FETS1D52B4ULRH

mp = MATS1D5Bond( Ef = 1., Af = 1.,
                    Em = 1., Am = 1.,
                    bond_fn = MFnLineArray(xdata = [0.,1.],
                                           ydata = [0.,1.]))  
                                               
fets_eval = FETS1D52B4ULRH(mats_eval = mp ) 

discr = ( 6, 1 )
# Discretization
fe_domain1 = FEGrid( coord_min = (0,0,0),
                           coord_max = (1.,.1,0.), 
                           shape   = discr,
                           n_nodal_dofs = fets_eval.n_nodal_dofs,
                           dof_r = fets_eval.dof_r,
                           geo_r = fets_eval.geo_r )

fe_domain2 = FEGrid( coord_min = (0.,.1,0),
                           coord_max = (1.,.2,0.), 
                           shape   = discr,
                           n_nodal_dofs = fets_eval.n_nodal_dofs,
                           dof_r = fets_eval.dof_r,
                           geo_r = fets_eval.geo_r )


fe_domain  = FEDomainList( subdomains = [ fe_domain1, fe_domain2] )

# Tseval for a discretized line domain
ts_eval1 = DOTSEval( fets_eval = fets_eval )
ts_eval2 = DOTSEval( fets_eval = fets_eval )
#ts_eval3 = DOTSEval( fets_eval = fets_eval )

ts_eval = DOTSListEval( dots_list = [ts_eval1, ts_eval2] )

ts = TS( tse = ts_eval,
         dof_resultants = True,
         sdomain = fe_domain,
         bcond_list =  [ BCDofGroup( var='u', value = 5.e-3, dims = [0],
                                     get_dof_method = fe_domain1.get_bottom_dofs ),
                         BCDofGroup( var='u', value = 0., dims = [0],
                                     get_dof_method = fe_domain2.get_bottom_dofs,
                                     get_link_dof_method = fe_domain1.get_top_dofs,
                                     link_coeffs = [1.] ),
                         BCDofGroup( var='u', value = 0., dims = [0],
                                     get_dof_method = fe_domain2.get_top_dofs )],                                                                        
         rtrace_list =  [ 
                         RTraceGraph(name = 'Fi,right over u_right (iteration)' ,
                               var_y = 'F_int', idx_y = 0,
                               var_x = 'U_k', idx_x = 0),
#                         RTraceDomainListField( name = 'Displacement', 
#                                                 var = 'u', 
#                                                 idx = 0 )
#                          RTraceDomainListField(name = 'Damage' ,
#                          var = 'omega', idx = 0,
#                          record_on = 'update',
#                          warp = False),
#                              RTraceDomainListField(name = 'Stress' ,
#                              var = 'sig_app', idx = 0,
#                              record_on = 'update',
#                              warp = False),
#                              RTraceDomainListField(name = 'Strain' ,
#                              var = 'eps_app', idx = 0,
#                              record_on = 'update',
#                              warp = False), 
                          ]             
        )

# Add the time-loop control
#global tloop
tloop = TLoop( tstepper = ts, tolerance = 1e-5,
                    tline  = TLine( min = 0.0,  step = 0.5, max = .5 ))

import cProfile
cProfile.run('tloop.eval()', 'tloop_prof' )

import pstats
p = pstats.Stats('tloop_prof')
p.strip_dirs()
print 'cumulative'
p.sort_stats('cumulative').print_stats(20)
print 'time'
p.sort_stats('time').print_stats(20)    
#tloop.eval()

from ibvpy.plugins.ibvpy_app import IBVPyApp
app = IBVPyApp( ibv_resource = tloop )
app.main()      
