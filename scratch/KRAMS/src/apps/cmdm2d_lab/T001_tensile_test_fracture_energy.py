'''
Example of a tensile test using the mats3d
'''

from ibvpy.api import \
    RTraceGraph, TLine
from ibvpy.mats.mats2D.mats2D_explore import \
    MATS2DExplore

from ibvpy.mats.mats2D.mats2D_cmdm.mats2D_cmdm import\
    MATS2DMicroplaneDamage, PhiFnStrainSoftening

from ibvpy.core.rtrace_eval import \
    RTraceEval

from time import time

from os.path import join

mats_eval = MATS2DMicroplaneDamage( E = 34000, nu = 0.25,
                                      model_version = 'stiffness',
                                      phi_fn = PhiFnStrainSoftening(
                                                                G_f = 0.001117,
                                                                f_t = 2.8968,
                                                                md = 0.0,
                                                                h = 1.0 ) )

mats2D_explore = \
    MATS2DExplore( mats2D_eval = mats_eval )
    
mats2D_explore.tloop.tstepper.rtrace_list =  \
                   [ RTraceGraph(name = 'stress - strain',
                                 var_x = 'eps_app', idx_x = 0,
                                 var_y = 'sig_app', idx_y = 0,
                                 record_on = 'update' ),
                     RTraceGraph(name = 'strain 0 - strain 1',
                                 var_x = 'eps_app', idx_x = 0,
                                 var_y = 'eps_app', idx_y = 1,
                                 record_on = 'update' ),
                     RTraceGraph(name = 'fracture energy time',
                                 var_x = 'time', idx_x = 0,
                                 var_y = 'fracture_energy', idx_y = 0,
                                 update_on = 'update' ) ]


#---------------------------
# calculation:
#---------------------------
tmax    = 0.01  #[m]
n_steps = 5
mats2D_explore.tloop.tline = TLine( min = 0.0,  step=tmax/n_steps, max = tmax )

mats2D_explore.tloop.eval()

from ibvpy.plugins.ibvpy_app import IBVPyApp
ibvpy_app = IBVPyApp( ibv_resource = mats2D_explore )
ibvpy_app.main()

