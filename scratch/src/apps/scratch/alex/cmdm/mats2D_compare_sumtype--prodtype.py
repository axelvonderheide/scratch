
from ibvpy.core.rtrace import \
    RTraceGraph,RTraceArraySnapshot

from mathkit.mfn.mfn_line.mfn_line import MFn1DDataGrid

if __name__ == '__main__':
    #--------------------------------------------------------------------------------
    # Example using the mats2d_explore 
    #--------------------------------------------------------------------------------
    from ibvpy.mats.mats2D.mats2D_explore import MATS2DExplore
    from ibvpy.mats.mats2D.mats2D_rtrace_cylinder import MATS2DRTraceCylinder
    from mats2D_cmdm_rtrace_Gf_mac import MATS2DMicroplaneDamageTraceGfmac
    from mats2D_cmdm_rtrace_Gf_mic import MATS2DMicroplaneDamageTraceGfmic
    from mats2D_cmdm import MA2DMicroplaneDamage
    
    mats2D_explore = \
        MATS2DExplore( mats2D_eval = MA2DMicroplaneDamage( elastic_debug = False ),
                       rtrace_list = [ RTraceGraph(name = 'strain 0 - stress 0',
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
                                       RTraceGraph(name = 'time - sig_norm',
                                                   var_x = 'time', idx_x = 0,
                                                   var_y = 'sig_norm', idx_y = 0,
                                                   record_on = 'update' ),
                                       RTraceGraph(name = 'time - phi_pdc',
                                                   var_x = 'time', idx_x = 0,
                                                   var_y = 'phi_pdc', idx_y = 0,
                                                   record_on = 'update' ),
                                       RTraceGraph(name = 'time - microplane damage',
                                                   var_x = 'time', idx_x = 0,
                                                   var_y = 'microplane_damage', idx_y = 0,
                                                   record_on = 'update' ),
                                       RTraceGraph(name = 'e_equiv - s_equiv',
                                                   var_x = 'equiv_microstrains', idx_x = 0,
                                                   var_y = 'equiv_microstresses', idx_y = 0,
                                                   record_on = 'update' ),

                                       # microplane fracture energy:
                                       MATS2DMicroplaneDamageTraceGfmic(name = 'time - G_f_micro',
                                                                     var_x = 'equiv_microstrains', idx_x = 0,
                                                                     var_y = 'equiv_microstresses', idx_y = 0,
                                                                     record_on = 'update' ),
                                       
                                       # macroscopic fracture energy:
                                       MATS2DMicroplaneDamageTraceGfmac(name = 'time - G_f_macro',
                                                                          var_x = 'eps_app', idx_x = 0,
                                                                          var_y = 'sig_app', idx_y = 0,
                                                                          record_on = 'update' ),


                                       MATS2DRTraceCylinder(name = 'Laterne',
                                                            var_axis    = 'time', idx_axis = 0,
                                                            var_surface = 'microplane_damage',
                                                            record_on = 'update' ),
                                       RTraceArraySnapshot(name = 'fracture energy contributions',
                                                           var = 'fracture_energy_list',
                                                           record_on = 'update' ),
                                       RTraceArraySnapshot(name = 'microplane damage',
                                                           var = 'microplane_damage',
                                                           record_on = 'update' ),
                                     ]
                       )

    from numpy import copy, arange
    
    mats2D_explore.mats2D_eval.model_version = 'compliance'

    # compare sum-type and product-type symmetrization for the 
    # compliance version for uniaxial tension

    ### 1 ###
    mats2D_explore.mats2D_eval.symmetrization = 'sum-type'
    mats2D_explore.tloop.eval()
    mp_se_ee = mats2D_explore.rtrace_mngr['e_equiv - s_equiv']
    mp_se_ee.redraw()
    e_equiv1 = copy( mp_se_ee.trace.xdata )
    s_equiv1 = copy( mp_se_ee.trace.ydata )


#    trace1 = MFn1DDataGrid( xdata = copy( mp_se_ee.trace.xdata ),
#                            ydata = copy( mp_se_ee.trace.ydata ) )


    ### 2 ###
    mats2D_explore.mats2D_eval.symmetrization = 'product-type'
    mats2D_explore.tloop.eval()
    mp_se_ee = mats2D_explore.rtrace_mngr['e_equiv - s_equiv']
    mp_se_ee.redraw()
    e_equiv2 = copy( mp_se_ee.trace.xdata )
    s_equiv2 = copy( mp_se_ee.trace.ydata )


#    ### 3 ###
#    mats2D_explore.mats2D_eval.symmetrization = 'product-type'
#    mats2D_explore.tloop.eval()
#    mp_se_ee = mats2D_explore.rtrace_mngr['e_equiv_projection - s_equiv']
#    mp_se_ee.redraw()
#    e_equiv_projection = copy( mp_se_ee.trace.xdata )
#    s_equiv3 = copy( mp_se_ee.trace.ydata )


    ### 2-1 ###
    xdata = arange( len( e_equiv1 ), dtype = float )
    diff_trace = MFn1DDataGrid( xdata = xdata,
                                ydata = e_equiv2 - e_equiv1 )
    diff_trace.configure_traits()


    ### 2-1 ###
    diff_trace = MFn1DDataGrid( xdata = xdata,
                                ydata = s_equiv2 - s_equiv1 )
    diff_trace.configure_traits()
    
#    ### 3-1 ###
#    diff_trace = MFn1DDataGrid( xdata = xdata,
#                                ydata = e_equiv_projection - e_equiv1 )
#    diff_trace.configure_traits()
    
    
#    from ibvpy.plugins.ibvpy_app import IBVPyApp
#    ibvpy_app = IBVPyApp( ibv_resource = mats2D_explore )
#    ibvpy_app.main()




