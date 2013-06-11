
# Example - uniaxial tension with weakened cross section using Masar's isotropic damage model 

from ibvpy.api import \
    TStepper as TS, RTraceGraph, RTraceDomainListField, TLoop, \
    TLine, BCDof, IBVPSolve as IS, DOTSEval
from ibvpy.mats.mats2D.mats2D_elastic.mats2D_elastic import MATS2DElastic
from ibvpy.api import BCDofGroup
from ibvpy.fets.fets2D.fets2D4q import FETS2D4Q
from ibvpy.fets.fets2D.fets2D9q import FETS2D9Q

from ibvpy.mats.mats2D.mats2D_elastic.mats2D_elastic import MATS2DElastic

from ibvpy.mats.mats2D.mats2D_sdamage.mats2D_sdamage import MATS2DScalarDamage
from ibvpy.mats.mats2D.mats2D_sdamage.strain_norm2d import Energy, Euclidean, Mises, Rankine, Mazars

from ibvpy.mats.mats2D.mats2D_cmdm.mats2D_cmdm import MATS2DMicroplaneDamage


from ibvpy.mats.mats_proxy import MATSProxy
from mathkit.geo.geo_ndgrid import GeoNDGrid
from mathkit.mfn.mfn_ndgrid.mfn_ndgrid import MFnNDGrid, GridPoint
from math import pi
from numpy import sin, cos, c_, arange, hstack

from time import time

#mp = MATS2DScalarDamage(E = 30000.,
#                        nu = 0.2,
#                        epsilon_0 = 0.005,
#                        epsilon_f = .5e-2,
#                        strain_norm_type = 'Mazars')
#
#mp = MATS2DMicroplaneDamage(E = 30000.,
#                          nu = 0.2,
#                          symmetrization = 'sum-type',
##                          n_mp = 7,
#                          model_version = 'compliance')

mp = MATS2DElastic(E = 30000.,
                   nu = 0.2,
                   stress_state = 'plane_stress')

#fets_eval = FETS2D4Q(mats_eval = mp) 
fets_eval = FETS2D9Q(mats_eval = mp) 

from ibvpy.mesh.fe_grid import FEGrid

# Geometric transformation
#
def arch_2d( points ):
    x = points[:,0]
    y = points[:,1]
    
    # outer radius of the arc
    R = 1.395
    # opening angle of the arc (=segment of a circle)
    beta = 100.18 * pi/180
    #thickness at the top of the arc
    D_top = 0.03
    #thickness at the bottom of the arc
    D_bottom = 0.04
    # thickness change between top and bottom
    D_delta = D_bottom - D_top
    # angle between horizontal line and starting point of the arc
    alpha = (pi - beta )/2

    # derive the length and heights of the carthesian grid
    L = x[-1] - x[0]
    H = y[-1] - y[0]

    # --------------------------------------------------
    # transformation from carthesian coordinates (x,y) 
    # to curvilinear (polar) coordinates (phi,r):
    # --------------------------------------------------
   
    # angle starts at the horizontal line turning right:
    phi = alpha + x / L * beta
 
    # variable thickness:
    # if phi <= (pi/2):
#    # (linear transition assumed)
#    D_var = D_bottom - (phi - alpha)/(beta/2.) * D_delta    
#    # if phi > (pi/2):
#    bool_phi = phi > pi/2
#    D_var[bool_phi] = D_top + (phi[bool_phi] - pi/2.)/(beta/2.) * D_delta    
    # (quadratic transition assumed)
    D_var = D_delta*((phi-alpha)/(pi/2-alpha))**2 - 2.0*D_delta*((phi-alpha)/(pi/2-alpha)) + D_bottom    
    
    # radius: 
    r = R - (H - y) / H * D_var
    
    # carthesian coordinates transformed in the arc geometry
    x,y = - r * cos( phi ), r * sin( phi ) 
    return c_[ x,y ]

#def arch_2d( points ):
#    x = points[:,0]
#    y = points[:,1]
#    l = x[-1] - x[0]
#    R = 10.
#    phi = x / l * pi
#    r = R + y
#    x,y = - r * cos( phi ), r * sin( phi ) 
#    return c_[ x,y ]


# The 'coord_max'-coordinates of FEGrid
# are derived in the method arch2d and 
# automatically considered in the transformation
length = 1.
height = 1.

# 100 linear elements
shape = (4*20,1)

# Discretization
domain = FEGrid( coord_max = (length,height,0.), 
                 shape   = shape,
                 fets_eval = fets_eval,
                 geo_transform = arch_2d )
#
#right_dofs, right_dof_r = domain.get_right_dofs()
#right_dof = right_dofs[0,0]
#print 'right_dof',right_dof

top_middle_dofs, top_middle_dof_r = domain.get_top_middle_dofs()
# dof in global y-direction
top_middle_dof = top_middle_dofs[0,1]
print 'top_middle_dofs',top_middle_dofs
print 'top_middle_dof' ,top_middle_dof

top_dofs, top_dof_r = domain.get_top_dofs()
# dof in global y-direction
print 'fets_eval.n_e_dofs' , fets_eval.n_e_dofs

if fets_eval.n_e_dofs == 8:
    top_left_middle_dofs = top_dofs[shape[0]/4, :]
elif fets_eval.n_e_dofs == 18:
    top_left_middle_dofs = top_dofs[shape[0]/2, :]
print 'top_left_middle_dofs' , top_left_middle_dofs

beta = 100.18 * pi/180
alpha = (pi - beta )/2
phi_left_middle = alpha + beta/4


tstepper = TS( 
     sdomain = domain,
     bcond_list =  [ # constraints at supports 
                     BCDofGroup( var='u', value = 0., dims = [0,1],
                                 get_dof_method = domain.get_left_dofs ),
                     BCDofGroup( var='u', value = 0., dims = [0,1],
                                 get_dof_method = domain.get_right_dofs ),
                                 
#                     # loading at L/4 (asymmetric loading case)                     
#                     BCDof( var='u', value =  cos(phi_left_middle), dof = top_left_middle_dofs[0] ), 
#                     BCDof( var='u', value = -sin(phi_left_middle), dof = top_left_middle_dofs[1] ),
                                         
                     # loading at top-middle of the arc
                     BCDofGroup( var='u', value = -1., dims = [1],
                              get_dof_method = domain.get_top_middle_dofs )
                    ],

     rtrace_list =  [
                 RTraceGraph(name = 'Fi,right over u_right (iteration)' ,

#                             # (asymmetric loading case):
#                             var_y = 'F_int', idx_y = top_left_middle_dofs[0,1],
#                             var_x = 'U_k'  , idx_x = top_left_middle_dofs[0,1],

                             # (loading at top-middle):
                             var_y = 'F_int', idx_y = top_middle_dofs[0,1],
                             var_x = 'U_k'  , idx_x = top_middle_dofs[0,1],
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
tloop = TLoop( tstepper = tstepper, tolerance = 1e-4,
               tline  = TLine( min = 0.0, step = 0.0015, max = 0.0045))                   
t1 = time()
tloop.eval()
t2 = time()
print 'time: ', t2 - t1

# Put the whole thing into the simulation-framework to map the
# individual pieces of definition into the user interface.
#
from ibvpy.plugins.ibvpy_app import IBVPyApp
app = IBVPyApp( ibv_resource = tloop )
app.main()    
