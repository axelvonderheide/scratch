#-------------------------------------------------------------------------------
#
# Copyright (c) 2009, IMB, RWTH Aachen.
# All rights reserved.
#
# This software is provided without warranty under the terms of the BSD
# license included in simvisage/LICENSE.txt and may be redistributed only
# under the conditions described in the aforementioned license.  The license
# is also available online at http://www.simvisage.com/licenses/BSD.txt
#
# Thanks for using Simvisage open source!
#
# Created on May 29, 2009 by: rchx

from enthought.traits.api import Float
from ibvpy.mats.mats2D.mats2D_explore import MATS2DExplore


class MATS2DDamageFitter( MATS2DExplore ):

    step_size = Float( 0.002 )

    def go_one_step(self):
        
        current_time = self.tloop.tline.val
        tmax = current_time + self.tloop.tline.step
        self.tloop.tline.min = current_time
        self.tloop.tline.max = tmax
        self.tloop.eval()


if __name__ == '__main__':
    #--------------------------------------------------------------------------------
    # Example using the mats2d_explore 
    #--------------------------------------------------------------------------------
    from ibvpy.mats.mats2D.mats2D_explore import MATS2DExplore
    from ibvpy.mats.mats2D.mats2D_rtrace_cylinder import MATS2DRTraceCylinder
    from ibvpy.mats.mats2D.mats2D_cmdm.mats2D_cmdm_rtrace_Gf_mic import MATS2DMicroplaneDamageTraceGfmic,\
        MATS2DMicroplaneDamageTraceEtmic, MATS2DMicroplaneDamageTraceUtmic
    from ibvpy.mats.mats2D.mats2D_cmdm.mats2D_cmdm_rtrace_Gf_mac import MATS2DMicroplaneDamageTraceGfmac,\
        MATS2DMicroplaneDamageTraceEtmac, MATS2DMicroplaneDamageTraceUtmac
    from ibvpy.mats.mats2D.mats2D_cmdm.mats2D_cmdm import MATS2DMicroplaneDamage   
    from ibvpy.api import RTraceGraph, RTraceArraySnapshot 
    
    fitter = MATS2DDamageFitter \
            (          mats2D_eval = MATS2DMicroplaneDamage( elastic_debug = False ),
                       rtrace_list = [ 
                                       MATS2DRTraceCylinder(name = 'Laterne',
                                                            var_axis    = 'time', idx_axis = 0,
                                                            var_surface = 'microplane_damage',
                                                            update_on = 'update' ),
                                                            
                                       RTraceGraph(name = 'strain 0 - stress 0',
                                                   var_x = 'eps_app_v3d', idx_x = 0,
                                                   var_y = 'sig_app_v3d', idx_y = 0,
                                                   record_on = 'update' ),
                                       RTraceGraph(name = 'strain 0 - strain 1',
                                                   var_x = 'eps_app_v3d', idx_x = 0,
                                                   var_y = 'eps_app_v3d', idx_y = 1,
                                                   record_on = 'update' ),
                                       RTraceGraph(name = 'stress 0 - stress 1',
                                                   var_x = 'sig_app_v3d', idx_x = 0,
                                                   var_y = 'sig_app_v3d', idx_y = 1,
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
                                                           
#                                       # G_f_mic: microplane fracture energy:
#                                       MATS2DMicroplaneDamageTraceGfmic(name = 'G_f_mic_equiv',
#                                                                        var_x = 'e_equiv_arr', idx_x = 0,
#                                                                        var_y = 's_equiv_arr', idx_y = 0,
#                                                                        record_on = 'update' ),
#                                       MATS2DMicroplaneDamageTraceGfmic(name = 'G_f_mic_N',
#                                                                        var_x = 'e_N_arr', idx_x = 0,
#                                                                        var_y = 's_N_arr', idx_y = 0,
#                                                                        record_on = 'update' ),
#                                       MATS2DMicroplaneDamageTraceGfmic(name = 'G_f_mic_T',
#                                                                        var_x = 'e_T_arr', idx_x = 0,
#                                                                        var_y = 's_T_arr', idx_y = 0,
#                                                                        record_on = 'update' ),
#                                       # E_t_mic: microplane total energy                                                                       
#                                       MATS2DMicroplaneDamageTraceEtmic(name = 'E_t_mic_equiv',
#                                                                        var_x = 'e_equiv_arr', idx_x = 0,
#                                                                        var_y = 's_equiv_arr', idx_y = 0,
#                                                                        record_on = 'update' ),
#                                       MATS2DMicroplaneDamageTraceEtmic(name = 'E_t_mic_N',
#                                                                        var_x = 'e_N_arr', idx_x = 0,
#                                                                        var_y = 's_N_arr', idx_y = 0,
#                                                                        record_on = 'update' ),
#                                       MATS2DMicroplaneDamageTraceEtmic(name = 'E_t_mic_T',
#                                                                        var_x = 'e_T_arr', idx_x = 0,
#                                                                        var_y = 's_T_arr', idx_y = 0,
#                                                                        record_on = 'update' ),
#                                       # U_t_mic: microplane elastic energy                                                                        
#                                       MATS2DMicroplaneDamageTraceUtmic(name = 'U_t_mic_equiv',
#                                                                        var_x = 'e_equiv_arr', idx_x = 0,
#                                                                        var_y = 's_equiv_arr', idx_y = 0,
#                                                                        record_on = 'update' ),
#                                       MATS2DMicroplaneDamageTraceUtmic(name = 'U_t_mic_N',
#                                                                        var_x = 'e_N_arr', idx_x = 0,
#                                                                        var_y = 's_N_arr', idx_y = 0,
#                                                                        record_on = 'update' ),
#                                       MATS2DMicroplaneDamageTraceUtmic(name = 'U_t_mic_T',
#                                                                        var_x = 'e_T_arr', idx_x = 0,
#                                                                        var_y = 's_T_arr', idx_y = 0,
#                                                                        record_on = 'update' ),
                                
                                       # direction 11:                  
                                       # G_f_mac: macroscopic fracture energy:
                                       MATS2DMicroplaneDamageTraceGfmac(name = 'G_f_mac_11',

                                                                        var_x = 'eps_app', idx_x = 0,
                                                                        var_y = 'sig_app', idx_y = 0,
                                                                        update_on = 'update' ),
                                       # E_t_mac: macroscopic total energy:
                                       MATS2DMicroplaneDamageTraceEtmac(name = 'E_t_mac_11',
                                                                        var_x = 'eps_app', idx_x = 0,
                                                                        var_y = 'sig_app', idx_y = 0,
                                                                        update_on = 'update' ),
                                       # U_t_mac: macroscopic elastic energy:
                                       MATS2DMicroplaneDamageTraceUtmac(name = 'U_t_mac_11',
                                                                        var_x = 'eps_app', idx_x = 0,
                                                                        var_y = 'sig_app', idx_y = 0,
                                                                        update_on = 'update' ),

#                                       # direction 22:
#                                       # G_f_mac: macroscopic fracture energy:
#                                       MATS2DMicroplaneDamageTraceGfmac(name = 'G_f_mac_22',
#                                                                        var_x = 'eps_app', idx_x = 1,
#                                                                        var_y = 'sig_app', idx_y = 1,
#                                                                        record_on = 'update' ),
#                                       # E_t_mac: macroscopic total energy:
#                                       MATS2DMicroplaneDamageTraceEtmac(name = 'E_t_mac_22',
#                                                                        var_x = 'eps_app', idx_x = 1,
#                                                                        var_y = 'sig_app', idx_y = 1,
#                                                                        record_on = 'update' ),
#                                       # U_t_mac: macroscopic elastic energy:
#                                       MATS2DMicroplaneDamageTraceUtmac(name = 'U_t_mac_22',
#                                                                        var_x = 'eps_app', idx_x = 1,
#                                                                        var_y = 'sig_app', idx_y = 1,
#                                                                        record_on = 'update' ),
#                                       
#                                       # direction 12:                                     
#                                       # G_f_mac: macroscopic fracture energy:
#                                       MATS2DMicroplaneDamageTraceGfmac(name = 'G_f_mac_12',
#                                                                        var_x = 'eps_app', idx_x = 2,
#                                                                        var_y = 'sig_app', idx_y = 2,
#                                                                        record_on = 'update' ),
#                                       # E_t_mac: macroscopic total energy:
#                                       MATS2DMicroplaneDamageTraceEtmac(name = 'E_t_mac_12',
#                                                                        var_x = 'eps_app', idx_x = 2,
#                                                                        var_y = 'sig_app', idx_y = 2,
#                                                                        record_on = 'update' ),
#                                       # U_t_mac: macroscopic elastic energy:
#                                       MATS2DMicroplaneDamageTraceUtmac(name = 'U_t_mac_12',
#                                                                        var_x = 'eps_app', idx_x = 2,
#                                                                        var_y = 'sig_app', idx_y = 2,
#                                                                        update_on = 'update' ),

#                                       RTraceArraySnapshot(name = 'fracture energy contributions',
#                                                           var = 'fracture_energy_arr',

#                                                           update_on = 'update' ),

#                                     ### decoupled energy contributions for G_f
                                     RTraceGraph(name = 'time - G_f',
                                                  var_x = 'time', idx_x = 0,
                                                  var_y = 'fracture_energy', idx_y = 0,
                                                  record_on = 'update' ),
#                                     ###
                                       RTraceGraph(name = 'time - sig_norm',
                                                   var_x = 'time', idx_x = 0,
                                                   var_y = 'sig_norm', idx_y = 0,
                                                   record_on = 'update' ),
#                                       RTraceGraph(name = 'time - phi_pdc',
#                                                   var_x = 'time', idx_x = 0,
#                                                   var_y = 'phi_pdc', idx_y = 0,
#                                                   record_on = 'update' ),
#                                       # e_equiv_projection:            
#                                       RTraceGraph(name = 'e_equiv_projection - s_equiv',
#                                                   var_x = 'equiv_projection', idx_x = 0,
#                                                   var_y = 's_equiv', idx_y = 0,
#                                                   record_on = 'update' ),

                                     ]
                       )


        
    for i in range(40):
        fitter.go_one_step()

    from ibvpy.plugins.ibvpy_app import IBVPyApp
    ibvpy_app = IBVPyApp( ibv_resource = fitter )
    ibvpy_app.main()    