
# Example - uniaxial tension with weakened cross section using Masar's isotropic damage model 

from ibvpy.api import \
    TStepper as TS, MGridDomain, RTraceGraph, RTraceDomainField, TLoop, \
    TLine, BCDof, IBVPSolve as IS, DOTSEval
from ibvpy.mats.mats2D.mats2D_elastic.mats2D_elastic import MATS2DElastic
from ibvpy.api import BCDofGroup
from ibvpy.fets.fets2D.fets2D4q import FETS2D4Q
from ibvpy.mats.mats2D.mats2D_sdamage.mats2D_sdamage import MATS2DScalarDamage
from ibvpy.mats.mats2D.mats2D_sdamage.strain_norm2d import Energy, Euclidean, Mises, Rankine, Mazars
from ibvpy.mats.mats_proxy import MATSProxy
from mathkit.geo.geo_ndgrid import GeoNDGrid
from mathkit.mfn.mfn_ndgrid.mfn_ndgrid import MFnNDGrid, GridPoint
from math import pi
from numpy import sin, cos, c_

length = 1.
height = 5.

mp = MATS2DScalarDamage(E = 1.,
                        nu = 0.,
                        epsilon_0 = 0.005,
                        epsilon_f = .5e-2,
                        strain_norm_type = 'Mazars')

fets_eval = FETS2D4Q(mats_eval = mp) 
tseval  = DOTSEval( fets_eval = fets_eval )

from ibvpy.mesh.fe_grid import FEGrid
from ibvpy.mesh.fe_level_set_domain import FELevelSetDomain

# Discretization

#fe_domain = FEGrid( coord_min = ( -2.,-2., 0 ),
#                          coord_max = (  2., 2., 0 ),
#                          shape = (5,5),
#                          n_nodal_dofs = 2,
#                          geo_r = [[-1,-1],
#                                        [-1, 1],
#                                        [ 1, 1],
#                                        [ 1,-1]],
#                          dof_r = [[-1,-1],
#                                        [-1, 1],
#                                        [ 1, 1],
#                                        [ 1,-1]] )
#    
#ls_domain = FELevelSetDomain( source_domain = fe_domain,
#                              ls_function = lambda x,y: x**2 + y**2 - 1. )
    

fe_domain = FEGrid( coord_min = ( -2.,-2., 0 ),
                          coord_max = (  2., 2., 0 ),
                          shape = (2,2),
                          n_nodal_dofs = 2,
                          geo_r = [[-1,-1],
                                        [-1, 1],
                                        [ 1, 1],
                                        [ 1,-1]],
                          dof_r = [[-1,-1],
                                        [-1, 1],
                                        [ 1, 1],
                                        [ 1,-1]] )
    
ls_domain = FELevelSetDomain( source_domain = fe_domain,
                              ls_function = lambda x,y: x )
    
#
tstepper = TS( tse = tseval, sdomain = ls_domain )

# Add the time-loop control
tloop = TLoop( tstepper = tstepper,
               DT = .00025,
               tline  = TLine( min = 0.0,  max = .0015))                   

K, F_int, F_ext = tloop.get_initial_state()
print K

# Put the whole thing into the simulation-framework to map the
# individual pieces of definition into the user interface.
#
from ibvpy.plugins.ibvpy_app import IBVPyApp
app = IBVPyApp( ibv_resource = tstepper )
app.main()    
