
# Example - uniaxial tension with weakened cross section using Masar's isotropic damage model 

from ibvpy.api import \
    TStepper as TS, RTraceGraph, RTraceDomainListField, TLoop, \
    TLine, BCDof, IBVPSolve as IS, DOTSEval
from ibvpy.mats.mats2D.mats2D_elastic.mats2D_elastic import MATS2DElastic
from ibvpy.api import BCDofGroup
from ibvpy.fets.fets2D.fets2D4q import FETS2D4Q
from ibvpy.mats.mats2D.mats2D_sdamage.mats2D_sdamage import \
    MATS2DScalarDamage
from ibvpy.mats.mats2D.mats2D_cmdm.mats2D_cmdm import \
    MATS2DMicroplaneDamage
from mathkit.geo.geo_ndgrid import GeoNDGrid
from mathkit.mfn.mfn_ndgrid.mfn_ndgrid import MFnNDGrid, GridPoint
from math import pi
from numpy import sin, cos, c_

length = 1.
height = 5.

mp = MATS2DMicroplaneDamage()

fets_eval = FETS2D4Q(mats_eval = mp) 

from ibvpy.mesh.fe_grid import FEGrid

# Geometric transformation
#
def arch_2d( points ):
    x = points[:,0]
    y = points[:,1]
    l = x[-1] - x[0]
    R = 10.
    phi = x / l * pi
    r = R + y
    x,y = - r * cos( phi ), r * sin( phi ) 
    return c_[ x,y ]

# Discretization
domain = FEGrid( coord_max = (length,height,0.), 
                       shape   = (10,5),
#                       shape   = (5,2),
                       fets_eval = fets_eval,
                       geo_transform = arch_2d )

#
right_dofs, right_dof_r = domain.get_right_dofs()
right_dof = right_dofs[0,0]
print 'right_dof',right_dof                         
tstepper = TS( sdomain = domain,
     bcond_list =  [ BCDofGroup( var='u', value = 0., dims = [0,1],
                              get_dof_method = domain.get_left_dofs ),
                     BCDofGroup( var='u', value = 0.01125, dims = [0],
                              get_dof_method = domain.get_right_dofs ) ],
     rtrace_list =  [
                 RTraceGraph(name = 'Fi,right over u_right (iteration)' ,
                             var_y = 'F_int', idx_y = right_dof,
                             var_x = 'U_k', idx_x = right_dof,
                             record_on = 'update'),
              RTraceDomainListField(name = 'Displacement' ,
              var = 'u', idx = 0,
              record_on = 'update',
              warp = True),
              RTraceDomainListField(name = 'Fracture Energy' ,
              var = 'fracture_energy', idx = 0,
              record_on = 'update',
              warp = True),
              RTraceDomainListField(name = 'Strain' ,
              var = 'eps_app', idx = 0,
              record_on = 'update',
              warp = False),
              RTraceDomainListField(name = 'Stress' ,
              var = 'sig_app', idx = 0,
              record_on = 'update',
              warp = False),                             
            ]             
        )

# Add the time-loop control
tloop = TLoop( tstepper = tstepper,
               tolerance = 1e-3,
               KMAX = 100,
               tline  = TLine( min = 0.0,  step = 0.2, max = 1.01))                   

tloop.eval()

# Put the whole thing into the simulation-framework to map the
# individual pieces of definition into the user interface.
#
from ibvpy.plugins.ibvpy_app import IBVPyApp
app = IBVPyApp( ibv_resource = tloop )
app.main()    
