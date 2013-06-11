
from ibvpy.core.rtrace import \
    RTraceGraph,RTraceArraySnapshot

from mathkit.mfn.mfn_line.mfn_line import MFn1DDataGrid

if __name__ == '__main__':
    #--------------------------------------------------------------------------------
    # Example using the mats2d_explore 
    #--------------------------------------------------------------------------------
    from ibvpy.mats.mats2D.mats2D_explore import MATS2DExplore
    from ibvpy.mats.mats2D.mats2D_rtrace_cylinder import MATS2DRTraceCylinder
    from mats2D_cmdm_rtrace_Gf_mic import MATS2DMicroplaneDamageTraceGfmic,\
        MATS2DMicroplaneDamageTraceEtmic, MATS2DMicroplaneDamageTraceUtmic
    from mats2D_cmdm_rtrace_Gf_mac import MATS2DMicroplaneDamageTraceGfmac,\
        MATS2DMicroplaneDamageTraceEtmac, MATS2DMicroplaneDamageTraceUtmac
    from mats2D_cmdm import MA2DMicroplaneDamage
    
    mats2D_explore = \
        MATS2DExplore( mats2D_eval = MA2DMicroplaneDamage( elastic_debug = False ),
                       rtrace_list = [ 
                                                                            
                                       # G_f_mic: microplane fracture energy:
                                       MATS2DMicroplaneDamageTraceGfmic(name = 'G_f_mic_equiv',
                                                                        var_x = 'e_equiv_arr', idx_x = 0,
                                                                        var_y = 's_equiv_arr', idx_y = 0,
                                                                        record_on = 'update' ),
                                       MATS2DMicroplaneDamageTraceGfmic(name = 'G_f_mic_N',
                                                                        var_x = 'e_N_arr', idx_x = 0,
                                                                        var_y = 's_N_arr', idx_y = 0,
                                                                        record_on = 'update' ),
                                       MATS2DMicroplaneDamageTraceGfmic(name = 'G_f_mic_T',
                                                                        var_x = 'e_T_arr', idx_x = 0,
                                                                        var_y = 's_T_arr', idx_y = 0,
                                                                        record_on = 'update' ),
                                       # E_t_mic: microplane total energy                                                                       
                                       MATS2DMicroplaneDamageTraceEtmic(name = 'E_t_mic_equiv',
                                                                        var_x = 'e_equiv_arr', idx_x = 0,
                                                                        var_y = 's_equiv_arr', idx_y = 0,
                                                                        record_on = 'update' ),
                                       MATS2DMicroplaneDamageTraceEtmic(name = 'E_t_mic_N',
                                                                        var_x = 'e_N_arr', idx_x = 0,
                                                                        var_y = 's_N_arr', idx_y = 0,
                                                                        record_on = 'update' ),
                                       MATS2DMicroplaneDamageTraceEtmic(name = 'E_t_mic_T',
                                                                        var_x = 'e_T_arr', idx_x = 0,
                                                                        var_y = 's_T_arr', idx_y = 0,
                                                                        record_on = 'update' ),
                                       # U_t_mic: microplane elastic energy                                                                        
                                       MATS2DMicroplaneDamageTraceUtmic(name = 'U_t_mic_equiv',
                                                                        var_x = 'e_equiv_arr', idx_x = 0,
                                                                        var_y = 's_equiv_arr', idx_y = 0,
                                                                        record_on = 'update' ),
                                       MATS2DMicroplaneDamageTraceUtmic(name = 'U_t_mic_N',
                                                                        var_x = 'e_N_arr', idx_x = 0,
                                                                        var_y = 's_N_arr', idx_y = 0,
                                                                        record_on = 'update' ),
                                       MATS2DMicroplaneDamageTraceUtmic(name = 'U_t_mic_T',
                                                                        var_x = 'e_T_arr', idx_x = 0,
                                                                        var_y = 's_T_arr', idx_y = 0,
                                                                        record_on = 'update' ),
  
                                   
                                       # direction 11:                  
                                       # G_f_mac: macroscopic fracture energy:
                                       MATS2DMicroplaneDamageTraceGfmac(name = 'G_f_mac_11',
                                                                        var_x = 'eps_app', idx_x = 0,
                                                                        var_y = 'sig_app', idx_y = 0,
                                                                        record_on = 'update' ),
                                       # E_t_mac: macroscopic total energy:
                                       MATS2DMicroplaneDamageTraceEtmac(name = 'E_t_mac_11',
                                                                        var_x = 'eps_app', idx_x = 0,
                                                                        var_y = 'sig_app', idx_y = 0,
                                                                        record_on = 'update' ),
                                       # U_t_mac: macroscopic elastic energy:
                                       MATS2DMicroplaneDamageTraceUtmac(name = 'U_t_mac_11',
                                                                        var_x = 'eps_app', idx_x = 0,
                                                                        var_y = 'sig_app', idx_y = 0,
                                                                        record_on = 'update' ),

                                       # direction 22:
                                       # G_f_mac: macroscopic fracture energy:
                                       MATS2DMicroplaneDamageTraceGfmac(name = 'G_f_mac_22',
                                                                        var_x = 'eps_app', idx_x = 1,
                                                                        var_y = 'sig_app', idx_y = 1,
                                                                        record_on = 'update' ),
                                       # E_t_mac: macroscopic total energy:
                                       MATS2DMicroplaneDamageTraceEtmac(name = 'E_t_mac_22',
                                                                        var_x = 'eps_app', idx_x = 1,
                                                                        var_y = 'sig_app', idx_y = 1,
                                                                        record_on = 'update' ),
                                       # U_t_mac: macroscopic elastic energy:
                                       MATS2DMicroplaneDamageTraceUtmac(name = 'U_t_mac_22',
                                                                        var_x = 'eps_app', idx_x = 1,
                                                                        var_y = 'sig_app', idx_y = 1,
                                                                        record_on = 'update' ),
                                       
                                       # direction 12:                                     
                                       # G_f_mac: macroscopic fracture energy:
                                       MATS2DMicroplaneDamageTraceGfmac(name = 'G_f_mac_12',
                                                                        var_x = 'eps_app', idx_x = 2,
                                                                        var_y = 'sig_app', idx_y = 2,
                                                                        record_on = 'update' ),
                                       # E_t_mac: macroscopic total energy:
                                       MATS2DMicroplaneDamageTraceEtmac(name = 'E_t_mac_12',
                                                                        var_x = 'eps_app', idx_x = 2,
                                                                        var_y = 'sig_app', idx_y = 2,
                                                                        record_on = 'update' ),
                                       # U_t_mac: macroscopic elastic energy:
                                       MATS2DMicroplaneDamageTraceUtmac(name = 'U_t_mac_12',
                                                                        var_x = 'eps_app', idx_x = 2,
                                                                        var_y = 'sig_app', idx_y = 2,
                                                                        record_on = 'update' ),


                           ]
                       )




    
    # ------------------------------------------------------------
    # settings for the calculation:    
    # ------------------------------------------------------------

    from numpy import copy


    # --- Settings for calculation: ---

    # loading specification:
    mats2D_explore.alpha_degree = 0.
#    mats2D_explore.alpha_degree = 30.
#    mats2D_explore.alpha_degree = 170.
    alpha_degree = mats2D_explore.alpha_degree

    mats2D_explore.mats2D_eval.model_version = 'compliance'
#    mats2D_explore.mats2D_eval.model_version = 'stiffness'
    mats2D_explore.mats2D_eval.symmetrization = 'product-type'
    mats2D_explore.mats2D_eval.double_constraint = False

    print '\n'


    
    # Combinations of variables to be traces:
    etype_list = ['G_f', 'E_t', 'U_t']
    vcomp_list = ['equiv', 'N', 'T']
    dir_list   = ['11','22','12']
    
    # ------------------------------------------------------------
    # MACROSCOPIC ENERGYS:     
    # ------------------------------------------------------------
    result_list_mac = {}
    result_list_mac['model_version'] = mats2D_explore.mats2D_eval.model_version 
    result_list_mac['symmetrization'] = mats2D_explore.mats2D_eval.symmetrization 
    result_list_mac['double_constraint'] = mats2D_explore.mats2D_eval.double_constraint 
    mats2D_explore.tloop.eval()

    for etype in etype_list:
        for dir in dir_list:
            
            rtrace_key = etype + '_mac_' + dir
            
            mp_se_ee = mats2D_explore.rtrace_mngr[ rtrace_key ]
            mp_se_ee.redraw()
    #        xdata = copy( mp_se_ee.trace.xdata )
            ydata = copy( mp_se_ee.trace.ydata )
            print 'ydata', ydata, '\n'
            result_list_mac[ rtrace_key ] = ydata[-1]
        
#    # --------------------------------------------------------------   
#    # CONSITENTLY DERIVED pairs of microplane strains and stresses:    
#    # --------------------------------------------------------------
#    # uses the same settings as has been used for the macroscopic values
#    result_list_mic = {}
#    result_list_mic['model_version'] = mats2D_explore.mats2D_eval.model_version 
#    result_list_mic['symmetrization'] = mats2D_explore.mats2D_eval.symmetrization 
#    result_list_mic['double_constraint'] = mats2D_explore.mats2D_eval.double_constraint 
#
#    for etype in etype_list:
#        for vcomp in vcomp_list:
#            
#            rtrace_key = etype + '_mic_' + vcomp
#            
#            mp_se_ee = mats2D_explore.rtrace_mngr[ rtrace_key ]
#            mp_se_ee.redraw()
##            xdata = copy( mp_se_ee.trace.xdata )
#            ydata = copy( mp_se_ee.trace.ydata )
#            # in this version of the tracer all values are the same:
#            result_list_mic[ rtrace_key ] = ydata[-1]
#
#    
#    # ------------------------------------------------------------
#    # DOUBLE CONSTRAINT based microplane strains and stresses:    
#    # ------------------------------------------------------------
#    # rerun with doublr constraint
#    mats2D_explore.mats2D_eval.double_constraint = True
#    result_list_mic_DC = {}
#    result_list_mic_DC['model_version'] = mats2D_explore.mats2D_eval.model_version 
#    result_list_mic_DC['symmetrization'] = mats2D_explore.mats2D_eval.symmetrization 
#    result_list_mic_DC['double_constraint'] = mats2D_explore.mats2D_eval.double_constraint 
#    mats2D_explore.tloop.eval()
#
#    for etype in etype_list:
#        for vcomp in vcomp_list:
#            
#            rtrace_key = etype + '_mic_' + vcomp
#            
#            mp_se_ee = mats2D_explore.rtrace_mngr[ rtrace_key ]
#            mp_se_ee.redraw()
##            xdata = copy( mp_se_ee.trace.xdata )
#            ydata = copy( mp_se_ee.trace.ydata )
#            result_list_mic_DC[ rtrace_key ] = ydata[-1]



    print '\n'
    print 'alpha_degree', alpha_degree
    print '\n'
    
    print '\n'
    print '--- MAC ---'
    # print the calculation settings:
    print 'symmetrization: ', result_list_mac['symmetrization']
    del result_list_mac['symmetrization']
    print 'model_version: ', result_list_mac['model_version']
    del result_list_mac['model_version']
    print 'double_constraint: ', result_list_mac['double_constraint']
    del result_list_mac['double_constraint']
    ## print the results:
    #for n in range(len( result_list_mac )):
    #    print result_list_mac.keys()[n] + " = %.4e" %(result_list_mac.values()[n])
    
#    
#    print '\n'
#    print '--- MIC ---'
#    # print the calculation settings:
#    print 'symmetrization: ', result_list_mic['symmetrization']
#    del result_list_mic['symmetrization']
#    print 'model_version: ', result_list_mic['model_version']
#    del result_list_mic['model_version']
#    print 'double_constraint: ', result_list_mic['double_constraint']
#    del result_list_mic['double_constraint']
#    ## print the results:
#    #for n in range(len( result_list_mic )):
#    #    print result_list_mic.keys()[n] + " = %.4e" %(result_list_mic.values()[n])
#    
#    
#    print '\n'
#    print '--- MIC_DC ---'
#    # print the calculation settings:
#    print 'symmetrization: ', result_list_mic_DC['symmetrization']
#    del result_list_mic_DC['symmetrization']
#    print 'model_version: ', result_list_mic_DC['model_version']
#    del result_list_mic_DC['model_version']
#    print 'double_constraint: ', result_list_mic_DC['double_constraint']
#    del result_list_mic_DC['double_constraint']
#    ## print the results:
#    #for n in range(len( result_list_mic_DC )):
#    #    print result_list_mic_DC.keys()[n] + " = %.4e" %(result_list_mic_DC.values()[n])
#    
#    
#    
##    # --- E_t ---
##    print '\n'
##    print '--- E_t ---'
##    
    print 'E_t_mac_11      ', result_list_mac['E_t_mac_11']
#    print 'E_t_mac_22      ', result_list_mac['E_t_mac_22']
#    print 'E_t_mac_12      ', result_list_mac['E_t_mac_12']
#    print 'E_t_mac_tot     ', result_list_mac['E_t_mac_11']+result_list_mac['E_t_mac_12']+result_list_mac['E_t_mac_22']
#    print '\n'
#    print 'E_t_mic_N       ', result_list_mic['E_t_mic_N']
#    print 'E_t_mic_T       ', result_list_mic['E_t_mic_T']
#    print 'E_t_mic_equiv   ', result_list_mic['E_t_mic_equiv']
#    print 'E_t_mic_tot     ', result_list_mic['E_t_mic_N'] + result_list_mic['E_t_mic_T']
#    print '\n'
#    print 'E_t_mic_N_DC    ', result_list_mic_DC['E_t_mic_N']
#    print 'E_t_mic_T_DC    ', result_list_mic_DC['E_t_mic_T']
#    print 'E_t_mic_equiv_DC', result_list_mic_DC['E_t_mic_equiv']
#    print 'E_t_mic_tot_DC  ', result_list_mic_DC['E_t_mic_N'] + result_list_mic_DC['E_t_mic_T']
#    
#    
#    # --- U_t ---
#    print '\n'
#    print '--- U_t ---'
#    
    print 'U_t_mac_11      ', result_list_mac['U_t_mac_11']
#    print 'U_t_mac_22      ', result_list_mac['U_t_mac_22']
#    print 'U_t_mac_12      ', result_list_mac['U_t_mac_12']
#    print 'U_t_mac_tot     ', result_list_mac['U_t_mac_11']+result_list_mac['U_t_mac_12']+result_list_mac['U_t_mac_22']
#    print '\n'
#    print 'U_t_mic_N       ', result_list_mic['U_t_mic_N']
#    print 'U_t_mic_T       ', result_list_mic['U_t_mic_T']
#    print 'U_t_mic_equiv   ', result_list_mic['U_t_mic_equiv']
#    print 'U_t_mic_tot     ', result_list_mic['U_t_mic_N'] + result_list_mic['U_t_mic_T']
#    print '\n'
#    print 'U_t_mic_N_DC    ', result_list_mic_DC['U_t_mic_N']
#    print 'U_t_mic_T_DC    ', result_list_mic_DC['U_t_mic_T']
#    print 'U_t_mic_equiv_DC', result_list_mic_DC['U_t_mic_equiv']
#    print 'U_t_mic_tot_DC  ', result_list_mic_DC['U_t_mic_N'] + result_list_mic_DC['U_t_mic_T']
#    
#    
#    
#    # --- G_f ---
#    print '\n'
#    print '--- G_f ---'
#    
    print 'G_f_mac_11      ', result_list_mac['G_f_mac_11']
#    print 'G_f_mac_22      ', result_list_mac['G_f_mac_22']
#    print 'G_f_mac_12      ', result_list_mac['G_f_mac_12']
#    print 'G_f_mac_tot     ', result_list_mac['G_f_mac_11']+result_list_mac['G_f_mac_12']+result_list_mac['G_f_mac_22']
#    print '\n'
#    print 'G_f_mic_N       ', result_list_mic['G_f_mic_N']
#    print 'G_f_mic_T       ', result_list_mic['G_f_mic_T']
#    print 'G_f_mic_equiv   ', result_list_mic['G_f_mic_equiv']
#    print 'G_f_mic_tot     ', result_list_mic['G_f_mic_N'] + result_list_mic['G_f_mic_T']
#    print '\n'
#    print 'G_f_mic_N_DC    ', result_list_mic_DC['G_f_mic_N']
#    print 'G_f_mic_T_DC    ', result_list_mic_DC['G_f_mic_T']
#    print 'G_f_mic_equiv_DC', result_list_mic_DC['G_f_mic_equiv']
#    print 'G_f_mic_tot_DC  ', result_list_mic_DC['G_f_mic_N'] + result_list_mic_DC['G_f_mic_T']
#    
#    # ---
