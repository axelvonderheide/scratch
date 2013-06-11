from ibvpy.core.api import \
    TStepper as TS, RTraceGraph, TLoop, \
    TLine, IBVPSolve as IS
    
from ibvpy.api import RTraceDomainField, BCDofGroup

from ibvpy.fets.fets2D.fets2D4q9u import  FETS2D4Q9U
from ibvpy.fets.fets2D.fets2D4q import FETS2D4Q
#from lib.fets.fets2D.fets2D4q8u import  FETS2D4Q8U
#from lib.fets.fets2D.fets2D9q import  FETS2D9Q
#from averaging import UniformDomainAveraging, LinearAF, QuarticAF
from ibvpy.api import DOTSEval
from ibvpy.mats.mats2D.mats2D_sdamage.mats2D_sdamage import MATS2DScalarDamage
from ibvpy.mats.mats2D.mats2D_sdamage.strain_norm2d import Energy, Euclidean, Mises, Rankine, Mazars
from ibvpy.mats.mats_proxy import MATSProxy
from mathkit.geo.geo_ndgrid import GeoNDGrid
from mathkit.mfn.mfn_ndgrid.mfn_ndgrid import MFnNDGrid, GridPoint
from mathkit.mfn.mfn_line.mfn_line import MFnLineArray

#fets_eval = FETS2D4Q9U(mats_eval = MATS2DScalarDamage(E = 34e3,
#                                               nu = 0.2,
#                                               epsilon_0 = 59e-6,
#                                               epsilon_f = 1.e-3,
#                                               #stiffness  = "algorithmic",
#                                               strain_norm = Mazars()))
#                                               #stiffness  = "algorithmic"))


#mp = MATSProxy( mats_eval_type = 'MATS2DScalarDamage')
mp = MATSProxy(mats_eval = MATS2DScalarDamage(E = 34.e3,
                               nu = 0.2,
                               epsilon_0 = 59.e-6,
                               epsilon_f = 3.2e-5,
                               #epsilon_f = 3.2e-1,
                               #stiffness  = "algorithmic",
                               strain_norm_type = 'Rankine'))
                                               #stiffness  = "algorithmic")))
mp.varpars['epsilon_0'].switch = 'varied'
mp.varpars['epsilon_0'].spatial_fn = MFnNDGrid( shape = (8,5,1),
                                    active_dims = ['x', 'y'],
                                    x_mins = GridPoint( x = 0., y = 0. ),
                                    x_maxs = GridPoint( x = 1., y = 0.5 ) )
mp.varpars['epsilon_0'].spatial_fn.set_values_in_box( 500., [0,0,0], [0.2,0.2,0.] )
mp.varpars['epsilon_0'].spatial_fn.set_values_in_box( 500., [0.8,0,0], [1.,0.2,0.] )
mp.varpars['epsilon_0'].spatial_fn.set_values_in_box( 50., [0.,0.46,0], [1.,0.5,0.] )


                                               
fets_eval = FETS2D4Q9U(mats_eval = mp)
#fets_eval = FETS2D4Q9U()
 
# Tseval for a discretized line domain
#
#tseval  = UniformDomainAveraging( fets_eval = fets_eval,
#                                  correction = True ,\
#                                  avg_function = QuarticAF(Radius = 0.2))
tseval  = DOTSEval( fets_eval = fets_eval)

# Discretization
#
from ibvpy.mesh.fe_grid import FEGrid
    
# Discretization
domain = FEGrid( coord_max = (1.,0.5,0.), 
                       shape   = (16, 8),
                       n_nodal_dofs = fets_eval.n_nodal_dofs,
                       dof_r = fets_eval.dof_r,
                       geo_r = fets_eval.geo_r )


mf = MFnLineArray( ydata = [0., 0.78, .79] )

# Put the tseval (time-stepper) into the spatial context of the
# discretization and specify the response tracers to evaluate there.
#

ts = TS( tse = tseval,
     sdomain = domain,
         bcond_list =  [BCDofGroup( var='u', value = 0., dims = [0,1],
                                  get_dof_method = domain.get_bottom_left_dofs ),
                        BCDofGroup( var='u', value = 0., dims = [1],
                                  get_dof_method = domain.get_bottom_right_dofs ),                                  
                        BCDofGroup( var='u', value = -1.2e-4, dims = [1],
                                    time_function = mf.get_value,
                                  get_dof_method = domain.get_top_middle_dofs ) ],
                    
     rtrace_list = [ 
#                        RTraceGraph(name = 'Fi,right over u_right (iteration)' ,
#                      var_y = 'F_int', idx_y = domain.get_bottom_left_dofs()[0,0],
#                      var_x = 'U_k', idx_x = domain.get_bottom_left_dofs()[0,0],
#                      record_on = 'update'),
#                      RTraceDomainField(name = 'Flux field' ,
#                      var = 'eps', idx = 0,
#                      record_on = 'update'),
              RTraceDomainField(name = 'Deformation' ,
              var = 'u', idx = 0,
              record_on = 'update',
              warp = True),
              RTraceDomainField(name = 'Damage' ,
              var = 'omega', idx = 0,
              record_on = 'update',
              warp = True),
              RTraceDomainField(name = 'Stress' ,
              var = 'sig_app', idx = 0,
              record_on = 'update',
              warp = False),
              RTraceDomainField(name = 'Strain' ,
              var = 'eps_app', idx = 0,
              record_on = 'update',
              warp = False),
              
#                      RTraceDomainField(name = 'N0' ,
#                                     var = 'N_mtx', idx = 2,
#                                     record_on = 'update')
#                      
         ]             
     )

# Add the time-loop control

tl = TLoop( tstepper = ts,
     tline  = TLine( min = 0.0, step = 1., max = 20. ),
     tolerance     =  1.e-6 )
if __name__ == '__main__':
    tl.eval()    
    # Put the whole stuff into the simulation-framework to map the
    # individual pieces of definition into the user interface.
    #
    from ibvpy.plugins.ibvpy_app import IBVPyApp
    app = IBVPyApp( ibv_resource = tl )
    app.main()