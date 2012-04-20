'''
Example of a tensile test using the mats2d_explore reducing
the 1d problem to an evaluation of a single material point
'''
from ibvpy.mats.mats2D.mats2D_explore import MATS2DExplore
from ibvpy.mats.mats2D.mats2D_rtrace_cylinder import MATS2DRTraceCylinder
from ibvpy.mats.mats2D.mats2D_cmdm.mats2D_cmdm_rtrace_Gf_mic import MATS2DMicroplaneDamageTraceGfmic,\
    MATS2DMicroplaneDamageTraceEtmic, MATS2DMicroplaneDamageTraceUtmic
from ibvpy.mats.mats2D.mats2D_cmdm.mats2D_cmdm_rtrace_Gf_mac import MATS2DMicroplaneDamageTraceGfmac,\
    MATS2DMicroplaneDamageTraceEtmac, MATS2DMicroplaneDamageTraceUtmac

from ibvpy.mats.mats2D.mats2D_cmdm.mats2D_cmdm import MATS2DMicroplaneDamage

from ibvpy.core.rtrace_eval import \
    RTraceEval

from ibvpy.api import \
    RTraceGraph,RTraceArraySnapshot

#from ibvpy.core.rtrace import \
#    RTraceGraph,RTraceArraySnapshot

from ibvpy.core.tloop import TLoop, TLine

from math import \
    pi as Pi, cos, sin, exp, sqrt as scalar_sqrt
    
from numpy import \
    sin, cos, c_, arange, hstack, array, loadtxt

from time import time

from os.path import join

mats2D_explore = \
    MATS2DExplore( mats2D_eval = MATS2DMicroplaneDamage( elastic_debug = False ),
                   rtrace_list = [ 
                                   MATS2DRTraceCylinder(name = 'Laterne',
                                                        var_axis    = 'time', idx_axis = 0,
                                                        var_surface = 'microplane_damage',
                                                        update_on = 'update' ),
                                                        
                                   RTraceGraph(name = 'strain 0 - stress 0',
                                               var_x = 'eps_app', idx_x = 0,
                                               var_y = 'sig_app', idx_y = 0,
                                               record_on = 'update' ),
                                   RTraceGraph(name = 'strain 0 - strain 1',
                                               var_x = 'eps_app', idx_x = 0,
                                               var_y = 'eps_app', idx_y = 1,
                                               record_on = 'update' ),
                                   RTraceGraph(name = 'stress 0 - stress 1',
                                               var_x = 'sig_app', idx_x = 0,
                                               var_y = 'sig_app', idx_y = 1,
                                               record_on = 'update' ),

                                   RTraceGraph(name = 'time - microplane damage',
                                               var_x = 'time', idx_x = 0,
                                               var_y = 'microplane_damage', idx_y = 0,
                                               update_on = 'update' ),

                                   # e_equiv, s_equiv
                                   RTraceGraph(name = 'e_equiv - s_equiv',
                                               var_x = 'e_equiv_arr', idx_x = 0,
                                               var_y = 's_equiv_arr', idx_y = 0,
                                               update_on = 'update' ),
                                
                                   # e_N, s_N:            
                                   RTraceGraph(name = 'e_N - s_N',
                                               var_x = 'e_N_arr', idx_x = 0,
                                               var_y = 's_N_arr', idx_y = 0,
                                               update_on = 'update' ),

                                   # e_T, s_T:            
                                   RTraceGraph(name = 'e_T - s_T',
                                               var_x = 'e_T_arr', idx_x = 0,
                                               var_y = 's_T_arr', idx_y = 0,
                                               update_on = 'update' ),

                                   RTraceArraySnapshot(name = 'equiv_projection',
                                                       var = 'equiv_projection',
                                                       record_on = 'update' ),
                                                   
                                   RTraceArraySnapshot(name = 'microplane damage',
                                                       var = 'microplane_damage',
                                                       record_on = 'update' ),

                                   RTraceArraySnapshot(name = 'e_equiv',
                                                       var = 'e_equiv_arr',
                                                       record_on = 'update' ),
                                   RTraceArraySnapshot(name = 's_equiv',
                                                       var = 's_equiv_arr',
                                                       record_on = 'update' ),
                                   RTraceArraySnapshot(name = 'e_N',
                                                       var = 'e_N_arr',
                                                       record_on = 'update' ),
                                   RTraceArraySnapshot(name = 's_N',
                                                       var = 's_N_arr',
                                                       record_on = 'update' ),
                                   RTraceArraySnapshot(name = 'e_T',
                                                       var = 'e_T_arr',
                                                       record_on = 'update' ),
                                   RTraceArraySnapshot(name = 's_T',
                                                       var = 's_T_arr',
                                                       record_on = 'update' ),
                                                       
                                   ###
                                   RTraceGraph(name = 'time - sig_norm',
                                               var_x = 'time', idx_x = 0,
                                               var_y = 'sig_norm', idx_y = 0,
                                               record_on = 'update' ),

                                     ]
                       )

        

# Elastic debug:    
mats2D_explore.mats2D_eval.elastic_debug = False
print 'elastic_debug: ', mats2D_explore.mats2D_eval.elastic_debug


#---------------------------
# Default settings:
#---------------------------

mats2D_explore.mats2D_eval.polar_fn_class = 'Isotropic polar function'
mats2D_explore.mats2D_eval.polar_fn.phi_fn_class = 'General'

# used data of phi_fn fitted from tensile test: 
fitted_phi_fn = loadtxt( join('mats_lab', 'fitted_phi_fn.out' ))
x = fitted_phi_fn[0]
y = fitted_phi_fn[1]
mats2D_explore.mats2D_eval.polar_fn.phi_fn.mfn.set( xdata = x, ydata = y )
mats2D_explore.mats2D_eval.polar_fn.phi_fn.mfn.data_changed = True

mats2D_explore.alpha_degree = 0.
print 'alpha_degree:  ', mats2D_explore.alpha_degree

mats2D_explore.mats2D_eval.stress_state = 'plane_stress'
print 'stress_state: ', mats2D_explore.mats2D_eval.stress_state

mats2D_explore.mats2D_eval.symmetrization = 'sum-type'
print 'symmetrization:', mats2D_explore.mats2D_eval.symmetrization

mats2D_explore.mats2D_eval.model_version = 'compliance'
print 'model_version: ', mats2D_explore.mats2D_eval.model_version


##---------------------------
## Anisotropic polar function - reinforced tensile test
##---------------------------
#mats2D_explore.mats2D_eval.polar_fn_class = 'Anisotropic polar function'
#
##    # Only available settings of 'AnisotropicPolarFn' for 'phi_fn_class':
##    mats2D_explore.mats2D_eval.polar_fn.phi_fn_class = 'TensionStiffening'
#
## default parameters of 'phi_fn':
#mats2D_explore.mats2D_eval.polar_fn.phi_fn.Epp    =  59e-6 # 0.5e-3     
#mats2D_explore.mats2D_eval.polar_fn.phi_fn.Efp    = 191e-6 # 4.0e-3     
#mats2D_explore.mats2D_eval.polar_fn.phi_fn.Elimit = 8.0e-3    
#mats2D_explore.mats2D_eval.polar_fn.phi_fn.Dfp = 0.4
#
## Polar function: parameter definition: 
#mats2D_explore.mats2D_eval.polar_fn.varied_params = ['Dfp']
#mats2D_explore.mats2D_eval.polar_fn.varpars['Dfp'].polar_fn.set( \
#                                                     phi_residual = 1.,
#                                                     phi_quasibrittle = 0.0,
#                                                     alpha = 0.,
#                                                     delta_trans = Pi/4. ,
#                                                     delta_alpha = Pi/8. ) 

#---------------------------
# calculation:
#---------------------------
tmax    = 0.01  #[m]
n_steps = 10
mats2D_explore.tloop.tline = TLine( min = 0.0,  step=tmax/n_steps, max = tmax )

#    mats2D_explore.mats2D_eval.configure_traits()
mats2D_explore.tloop.eval()
mats2D_explore.tloop.setup()

from ibvpy.plugins.ibvpy_app import IBVPyApp
ibvpy_app = IBVPyApp( ibv_resource = mats2D_explore )
ibvpy_app.main()

