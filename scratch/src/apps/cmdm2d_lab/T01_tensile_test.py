'''
Example of a tensile test using a one element discretization
'''
from ibvpy.api import \
    TStepper as TS, RTraceGraph, RTraceDomainListField, TLoop, \
    TLine, BCDof, IBVPSolve as IS, BCDofGroup

from ibvpy.fets.fets2D.fets2D4q import FETS2D4Q
from ibvpy.fets.fets2D.fets2D9q import FETS2D9Q

from ibvpy.mesh.fe_grid import FEGrid

from ibvpy.mats.mats2D.mats2D_elastic.mats2D_elastic import \
    MATS2DElastic
from ibvpy.mats.mats2D.mats2D_sdamage.mats2D_sdamage import \
    MATS2DScalarDamage
from ibvpy.mats.mats2D.mats2D_sdamage.strain_norm2d import \
    Energy, Euclidean, Mises, Rankine, Mazars
from ibvpy.mats.mats2D.mats2D_cmdm.mats2D_cmdm import \
    MATS2DMicroplaneDamage, PhiFnGeneral, PhiFnStrainSoftening   
from ibvpy.mats.mats_proxy import \
    MATSProxy
    
from mathkit.geo.geo_ndgrid import \
    GeoNDGrid
from mathkit.mfn.mfn_ndgrid.mfn_ndgrid import \
    MFnNDGrid, GridPoint

from numpy import sin, cos, c_, arange, hstack, array, loadtxt

from time import time
from os.path import join

from math import \
    pi as Pi, cos, sin, exp, sqrt as scalar_sqrt



#-------------------------
# calculation parameters:
#-------------------------

# geometry:
length = 0.45
height = 0.1

# one element discretization:
shape = (20,1)

# loading:
t_max   = 0.005
n_steps = 20
tline = TLine( min = 0.0,  step = t_max/n_steps, max = t_max)

element_type = 'linear'
#element_type = 'quadratic'


#-------------------------
# material model:
#-------------------------

read_fitted = True

if read_fitted:
    file_name = 'mats_lab/9-a_MAG-03-07.mats'
    file = open( file_name, 'r' )
    import pickle
    mfn = pickle.load( file )
    
    thickness = 0.03
    cmdm = MATS2DMicroplaneDamage( 
                                  E = 34000 * thickness,
                                  nu = 0.25,                                       
                                  n_mp = 30,
                                  elastic_debug = False,
                                  stress_state = 'plane_stress',
                                  symmetrization = 'sum-type',
                                  model_version = 'compliance',
                                  phi_fn = PhiFnGeneral( mfn = mfn ),
                                  )
    
else:
    cmdm = MATS2DMicroplaneDamage(E = 34000. * 0.04,
                            nu = 0.2,
                            symmetrization = 'sum-type',
                            model_version = 'compliance',
                            )
    cmdm.polar_fn_class = 'Isotropic polar function'
    cmdm.polar_fn.phi_fn_class = 'General'

    # used data of phi_fn fitted from tensile test: 
    fitted_phi_fn = loadtxt( join('mats_lab', 'fitted_phi_fn.out' ))
    x = fitted_phi_fn[0]
    y = fitted_phi_fn[1]
    cmdm.polar_fn.phi_fn.mfn.set( xdata = x, ydata = y )
    cmdm.polar_fn.phi_fn.mfn.data_changed = True
    
 
#-------------------------
# element type:
#-------------------------
if element_type == 'linear':
    fets_eval = FETS2D4Q( mats_eval = cmdm ) 
elif element_type == 'quadratic':
    fets_eval = FETS2D9Q( mats_eval = cmdm) 


#-------------------------
# Discretization
#-------------------------

domain = FEGrid( coord_max = (length,height,0.), 
                 shape   = shape,
                 fets_eval = fets_eval )

# OSOLET: Use individual domains instead for notched specimen. 
# For a notched beam deactivate the center element at the bottom of the beam:
# domain.deactivate( (shape[0]/2,0) )
# Note: For the center element an odd numer is necessary for n_elem[0].
#       For the center loading at top an even numkber is necessary.
       

#-------------------------
# ts: 
#-------------------------

#right_dofs = domain[0, 0, 0, 0].dofs[0,0,0]
#print 'right_dofs', right_dofs


right_dofs, right_dofs_points = domain.get_right_dofs()
print 'right_dofs', right_dofs

right_bottom_dofs, right_bottom_dofs_points = domain.get_bottom_right_dofs()
print 'right_bottom_dofs', right_bottom_dofs

tstepper = TS( dof_resultants = True,
               sdomain = domain,
     bcond_list =  [ # loading at the right hand side: 
                     BCDofGroup( var='u', value = 1.0, dims = [0],
                              get_dof_method = domain.get_right_dofs ),

                     # supports at the left hand side:
                     BCDofGroup( var='u', value = 0., dims = [0],
                              get_dof_method = domain.get_left_dofs ),
                     BCDofGroup( var='u', value = 0., dims = [0,1],
                              get_dof_method = domain.get_bottom_left_dofs ) ],
     rtrace_list = [
                     RTraceGraph(name = 'Fi at right over u at right (iteration)' ,
                                 var_y = 'F_int', idx_y = right_bottom_dofs[0,0],
                                 var_x = 'U_k'  , idx_x = right_bottom_dofs[0,0],
                                 record_on = 'update'),
                     
                     RTraceDomainListField(name = 'Displacement' ,
                     var = 'u', idx = 0,
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
tloop = TLoop( tstepper = tstepper, tolerance = 1e-4,
               tline  = tline )                   

tloop.eval()

# Put the whole thing into the simulation-framework to map the
# individual pieces of definition into the user interface.
#
from ibvpy.plugins.ibvpy_app import IBVPyApp
app = IBVPyApp( ibv_resource = tloop )
app.main()    

