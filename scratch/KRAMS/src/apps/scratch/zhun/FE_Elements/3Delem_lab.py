n_nod = 8
if n_nod == 27:
    dof_r =  \
                [[-1,-1,-1],
                 [ 1,-1,-1],
                 [ 1, 1,-1],
                 [-1, 1,-1],
                 [ 0,-1,-1],
                 [ 1, 0,-1],
                 [ 0, 1,-1],
                 [-1, 0,-1],
                 [-1,-1, 0],                   
                 [ 1,-1, 0],
                 [ 1, 1, 0],
                 [-1, 1, 0],
                 [-1,-1, 1],
                 [ 1,-1, 1],
                 [ 1, 1, 1],
                 [-1, 1, 1],
                 [ 0,-1, 1],
                 [ 1, 0, 1],
                 [ 0, 1, 1],
                 [-1, 0, 1],
                 [ 0, 0,-1],
                 [ 0,-1, 0],
                 [ 1, 0, 0],
                 [ 0, 1, 0],
                 [-1, 0, 0],
                 [ 0, 0, 1],
                 [ 0, 0, 0]]           
    geo_r =  \
         [[ -1., -1., -1.],
          [  1., -1., -1.],
          [ -1.,  1., -1.],
          [  1.,  1., -1.],
          [ -1., -1.,  1.],
          [  1., -1.,  1.],
          [ -1.,  1.,  1.],
          [  1.,  1.,  1.]]

elif n_nod == 20:
    geo_r = \
             [[ -1., -1., -1.],
              [  1., -1., -1.],
              [ -1.,  1., -1.],
              [  1.,  1., -1.],
              [ -1., -1.,  1.],
              [  1., -1.,  1.],
              [ -1.,  1.,  1.],
              [  1.,  1.,  1.]]

    dof_r = \
            [  [-1,-1,-1],
               [ 1,-1,-1],
               [ 1, 1,-1],
               [-1, 1,-1],
               [ 0,-1,-1],
               [ 1, 0,-1],
               [ 0, 1,-1],
               [-1, 0,-1],
               [-1,-1, 0],                   
               [ 1,-1, 0],
               [ 1, 1, 0],
               [-1, 1, 0],
               [-1,-1, 1],
               [ 1,-1, 1],
               [ 1, 1, 1],
               [-1, 1, 1],
               [ 0,-1, 1],
               [ 1, 0, 1],
               [ 0, 1, 1],
               [-1, 0, 1]]   

elif n_nod == 8:
    dof_r = \
             [[ -1., -1., -1.],
              [  1., -1., -1.],
              [ -1.,  1., -1.],
              [  1.,  1., -1.],
              [ -1., -1.,  1.],
              [  1., -1.,  1.],
              [ -1.,  1.,  1.],
              [  1.,  1.,  1.]]

    geo_r = \
             [[ -1., -1., -1.],
              [  1., -1., -1.],
              [ -1.,  1., -1.],
              [  1.,  1., -1.],
              [ -1., -1.,  1.],
              [  1., -1.,  1.],
              [ -1.,  1.,  1.],
              [  1.,  1.,  1.]]

from ibvpy.api import \
    TStepper as TS, MGridDomain, RTraceGraph, RTraceDomainField, TLoop, \
    TLine, BCDof, IBVPSolve as IS, DOTSEval
    
#from lib.mats.mats2D.mats_cmdm2D.mats_mdm2d import MACMDM
from ibvpy.mats.mats2D.mats2D_sdamage.mats2D_sdamage import MATS2DScalarDamage
from ibvpy.mats.mats2D.mats2D_sdamage.strain_norm2d import *
from ibvpy.mats.mats3D.mats3D_elastic.mats3D_elastic import MATS3DElastic
from ibvpy.fets.fets3D.fets3D8h27u import FETS3D8H27U
from ibvpy.fets.fets3D.fets3D8h20u import FETS3D8H20U
from ibvpy.fets.fets3D.fets3D8h import FETS3D8H

from ibvpy.mesh.mgrid_domain import MeshGridAdaptor
from ibvpy.mesh.fe_grid import FEGrid

from ibvpy.api import BCDofGroup
from mathkit.mfn.mfn_line.mfn_line import MFnLineArray
from ibvpy.mats.mats3D.mats3D_sdamage.mats3D_sdamage import MATS3DScalarDamage
from ibvpy.mats.mats3D.mats3D_sdamage.strain_norm3d import Euclidean,Energy,\
Mazars, Rankine, Mises

mats = MATS3DElastic(E=1.,nu=0.3)
#mats = MATS3DScalarDamage(E = 30000,
#                          nu = 0.2,
#                          epsilon_0 = 2.e-4,
#                          epsilon_f = 8.e-4,
#                          #strain_norm = Energy())
#                          #strain_norm = Euclidean())
#                          strain_norm = Mazars())

if n_nod == 27:
    fets_eval = FETS3D8H27U(mats_eval = mats) 
elif n_nod == 20:
    fets_eval = FETS3D8H20U(mats_eval = mats)  
elif n_nod == 8:
    fets_eval = FETS3D8H(mats_eval = mats)  

# Tseval for a discretized line domain
#
tseval  = DOTSEval( fets_eval = fets_eval )
 
# Discretization
#
domain = FEGrid( coord_max = (0.004,0.004,0.006), 
                       shape   = (3,3,3),
                       n_nodal_dofs = 3,
                       dof_r =  dof_r,           
                       geo_r =  geo_r)
                       
                       
# Put the tseval (time-stepper) into the spatial context of the
# discretization and specify the response tracers to evaluate there.
#

mf = MFnLineArray( ydata = [0.,0.0001,0.00015] )

ts = TS(
        tse = tseval,
        sdomain = domain,
         # conversion to list (square brackets) is only necessary for slicing of 
         # single dofs, e.g "get_left_dofs()[0,1]" which elsewise retuns an integer only
         bcond_list =  [ # constraint for all left dofs in x-direction: 
                     BCDofGroup(var='u', value = 0.,dims = [2],
                            get_dof_method = domain.get_back_dofs),\
                      # imposed displacement for all right dofs in y-direction:
                     BCDofGroup(var='u', value = 0.,dims = [1],
                            get_dof_method = domain.get_bottom_dofs),\
                     BCDofGroup(var='u', value = 0.,dims = [0],
                            get_dof_method = domain.get_left_dofs),\
#                     BCDofGroup(var='u', value = 0.,dims = [0,1,2],
#                            get_dof_method = domain.get_right_dofs),\
                     BCDofGroup(var='u', value = 1., dims = [2],
                             time_function = mf.get_value,
                            get_dof_method = domain.get_left_front_bottom_dof )]
                    ,
         rtrace_list = [ 
#                        RTraceGraph(name = 'Fi,right over u_right (iteration)' ,
#                                  var_y = 'F_int', idx_y = right_dof,
#                                  var_x = 'U_k', idx_x = right_dof),
#                        RTraceDomainField(name = 'Damage' ,
#                                          #position = 'int_pnts',
#                                          var = 'omega', idx = 0),
                        RTraceDomainField(name = 'Deformation x' ,
                                          #position = 'int_pnts',
                                          var = 'eps', idx = 0),
                        RTraceDomainField(name = 'Displacement' ,
                                          var = 'u', idx = 0),
                        RTraceDomainField(name = 'Stress' ,
                                          var = 'sig', idx = 0)

                    ]             
        )

# Add the time-loop control
#
tloop = TLoop( tstepper = ts,
         tline  = TLine( min = 0.0, step = 0.5, max = 2. ))
  
tloop.eval()    

# Put the whole stuff into the simulation-framework to map the
# individual pieces of definition into the user interface.
#
from ibvpy.plugins.ibvpy_app import IBVPyApp
app = IBVPyApp( ibv_resource = tloop )
app.main()
