
# Example - uniaxial tension with weakened cross section using Masar's isotropic damage model 

from ibvpy.api import \
    TStepper as TS, RTraceGraph, RTraceDomainListField, TLoop, \
    TLine, BCDof, IBVPSolve as IS
from ibvpy.mats.mats2D.mats2D_elastic.mats2D_elastic import MATS2DElastic
from ibvpy.api import BCDofGroup
from ibvpy.fets.fets2D.fets2D4q import FETS2D4Q
from ibvpy.mats.mats2D.mats2D_sdamage.mats2D_sdamage import MATS2DScalarDamage
from ibvpy.mats.mats2D.mats2D_sdamage.strain_norm2d import Energy, Euclidean, Mises, Rankine, Mazars
from ibvpy.mats.mats_proxy import MATSProxy
from mathkit.geo.geo_ndgrid import GeoNDGrid
from mathkit.mfn.mfn_ndgrid.mfn_ndgrid import MFnNDGrid, GridPoint

#fets_eval = FETS2D4Q9U(mats_eval = MATS2DScalarDamage(E = 34e3,
#                                               nu = 0.2,
#                                               epsilon_0 = 59e-6,
#                                               epsilon_f = 1.e-3,
#                                               #stiffness  = "algorithmic",
#                                               strain_norm = Mazars()))
#                                               #stiffness  = "algorithmic"))

length = 1.
height = 0.5

#mp = MATSProxy( mats_eval_type = 'MATS2DScalarDamage')
mp = MATSProxy( mats_eval = MATS2DScalarDamage(E = 1.,
                                               nu = 0.,
                                               epsilon_0 = 0.005,
                                               epsilon_f = .5e-2,
                                               #stiffness  = "algorithmic",
                                               strain_norm_type = 'Mazars'))
                                               #stiffness  = "algorithmic")))
mp.varpars['epsilon_0'].switch = 'varied'
mp.varpars['epsilon_0'].spatial_fn = MFnNDGrid( shape = (100,50,1),
                                    active_dims = ['x', 'y'],
                                    x_mins = GridPoint( x = 0., y = 0. ),
                                    x_maxs = GridPoint( x = length, y = height ) )
#mp.varpars['epsilon_0'].spatial_fn.set_values_in_box( 10., [0,0,0], [1.,0.5,0.] )
mp.varpars['epsilon_0'].spatial_fn.set_values_in_box( .1, [0.4,0,0], [0.5,0.5,0.] )
#mp.varpars['epsilon_0'].spatial_fn.set_values_in_box( 500., [0.8,0,0], [1.,0.2,0.] )
#mp.varpars['epsilon_0'].spatial_fn.set_values_in_box( 50., [0.,0.46,0], [1.,0.5,0.] )


fets_eval = FETS2D4Q(mats_eval = MATS2DElastic()) 
#fets_eval = FETS2D4Q(mats_eval = mp) 

from ibvpy.mesh.fe_grid import FEGrid

# Discretization
domain = FEGrid( coord_max = (length,height,0.), 
                       shape   = (10,5),
                       fets_eval = fets_eval )

#
right_dofs, right_dof_r = domain.get_right_dofs()
right_dof = right_dofs[0,0]
print 'right_dof',right_dof                         
tstepper = TS( 
     sdomain = domain,
     bcond_list =  [ BCDofGroup( var='u', value = 0., dims = [0,1],
                              get_dof_method = domain.get_left_dofs ),
                     BCDofGroup( var='u', value = 1., dims = [0],
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
#              RTraceDomainListField(name = 'Damage' ,
#              var = 'omega', idx = 0,
#              record_on = 'update',
#              warp = True),
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
               tline  = TLine( min = 0.0,  step = 0.00025, max = .0015))                   

tloop.eval()

# Put the whole thing into the simulation-framework to map the
# individual pieces of definition into the user interface.
#
from ibvpy.plugins.ibvpy_app import IBVPyApp
app = IBVPyApp( ibv_resource = tloop )
app.main()    
