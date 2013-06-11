
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
                                       # e_N, s_N:            
                                       RTraceGraph(name = 'MP0: e_N - s_N',
                                                   var_x = 'e_N_arr', idx_x = 0,
                                                   var_y = 's_N_arr', idx_y = 0,
                                                   record_on = 'update' ),
    
                                       # e_N, s_N:            
                                       RTraceGraph(name = 'MP1: e_N - s_N',
                                                   var_x = 'e_N_arr', idx_x = 1,
                                                   var_y = 's_N_arr', idx_y = 1,
                                                   record_on = 'update' ),
                                       
                                       # e_N, s_N:            
                                       RTraceGraph(name = 'MP2: e_N - s_N',
                                                   var_x = 'e_N_arr', idx_x = 2,
                                                   var_y = 's_N_arr', idx_y = 2,
                                                   record_on = 'update' ),
                           
                                       # e_N, s_N:            
                                       RTraceGraph(name = 'MP3: e_N - s_N',
                                                   var_x = 'e_N_arr', idx_x = 3,
                                                   var_y = 's_N_arr', idx_y = 3,
                                                   record_on = 'update' ),
                           
                           ]
                       )
    
    # ------------------------------------------------------------
    # settings for the calculation:    
    # ------------------------------------------------------------

    from numpy import copy, vstack
    from enthought.traits.api import Float
    
    # loading specification
    mats2D_explore.alpha_degree = 0.
    alpha_degree = mats2D_explore.alpha_degree
    mats2D_explore.mats2D_eval.model_version = 'stiffness'
    mats2D_explore.mats2D_eval.symmetrization = 'sum-type'


    
    ## no.1: 'e_N - s_N' - settings: ###
    mats2D_explore.mats2D_eval.double_constraint = False
    mversion_1    = mats2D_explore.mats2D_eval.model_version
    stype_1       = mats2D_explore.mats2D_eval.symmetrization
    dconstraint_1 = mats2D_explore.mats2D_eval.double_constraint
    mats2D_explore.tloop.eval()
    # MP1:
    mp_se_ee = mats2D_explore.rtrace_mngr['MP0: e_N - s_N']
    mp_se_ee.redraw()
    xdata_0 = copy( mp_se_ee.trace.xdata )
    ydata_0 = copy( mp_se_ee.trace.ydata )
    # MP2:
    mp_se_ee = mats2D_explore.rtrace_mngr['MP1: e_N - s_N']
    mp_se_ee.redraw()
    xdata_1 = copy( mp_se_ee.trace.xdata )
    ydata_1 = copy( mp_se_ee.trace.ydata )
    # MP3:
    mp_se_ee = mats2D_explore.rtrace_mngr['MP2: e_N - s_N']
    mp_se_ee.redraw()
    xdata_2 = copy( mp_se_ee.trace.xdata )
    ydata_2 = copy( mp_se_ee.trace.ydata )
    # MP4:
    mp_se_ee = mats2D_explore.rtrace_mngr['MP3: e_N - s_N']
    mp_se_ee.redraw()
    xdata_3 = copy( mp_se_ee.trace.xdata )
    ydata_3 = copy( mp_se_ee.trace.ydata )    
    # total load history:
#    e_N_no1 = vstack([xdata_0, xdata_1, xdata_2, xdata_3])
#    s_N_no1 = vstack([ydata_0, ydata_1, ydata_2, ydata_3])
    # last step:
    e_N_no1 = vstack([xdata_0[-1], xdata_1[-1], xdata_2[-1], xdata_3[-1]])
    s_N_no1 = vstack([ydata_0[-1], ydata_1[-1], ydata_2[-1], ydata_3[-1]])


    ## no2: 'e_N - s_N' - settings: ###
    mats2D_explore.mats2D_eval.double_constraint = True
#    mats2D_explore.mats2D_eval.model_version = 'stiffness'
    mversion_2    = mats2D_explore.mats2D_eval.model_version
    stype_2       = mats2D_explore.mats2D_eval.symmetrization
    dconstraint_2 = mats2D_explore.mats2D_eval.double_constraint
    mats2D_explore.tloop.eval()
    # MP1:
    mp_se_ee = mats2D_explore.rtrace_mngr['MP0: e_N - s_N']
    mp_se_ee.redraw()
    xdata_0 = copy( mp_se_ee.trace.xdata )
    ydata_0 = copy( mp_se_ee.trace.ydata )
    # MP2:
    mp_se_ee = mats2D_explore.rtrace_mngr['MP1: e_N - s_N']
    mp_se_ee.redraw()
    xdata_1 = copy( mp_se_ee.trace.xdata )
    ydata_1 = copy( mp_se_ee.trace.ydata )
    # MP3:
    mp_se_ee = mats2D_explore.rtrace_mngr['MP2: e_N - s_N']
    mp_se_ee.redraw()
    xdata_2 = copy( mp_se_ee.trace.xdata )
    ydata_2 = copy( mp_se_ee.trace.ydata )
    # MP4:
    mp_se_ee = mats2D_explore.rtrace_mngr['MP3: e_N - s_N']
    mp_se_ee.redraw()
    xdata_3 = copy( mp_se_ee.trace.xdata )
    ydata_3 = copy( mp_se_ee.trace.ydata )    
    # total load history:
#    e_N_no2 = vstack([xdata_0, xdata_1, xdata_2, xdata_3])
#    s_N_no2 = vstack([ydata_0, ydata_1, ydata_2, ydata_3])
    # last step:
    e_N_no2 = vstack([xdata_0[-1], xdata_1[-1], xdata_2[-1], xdata_3[-1]])
    s_N_no2 = vstack([ydata_0[-1], ydata_1[-1], ydata_2[-1], ydata_3[-1]])
    
    print '--- Settings for calculation: ---'
    print 'alpha_degree', alpha_degree
    print '\n'
    print 'model_version', mversion_1
    print 'symmetrization', stype_1
    print 'double_constraint', dconstraint_1
#    print 'e_N', e_N
    print '\n'
    print 'model_version', mversion_2
    print 'symmetrization', stype_2
    print 'double_constraint', dconstraint_2
#    print 'e_N_DC', e_N_DC
    print '\n'
    print 'e_N_no2 - e_N_no1', e_N_no2 - e_N_no1
    print 's_N_no2 - s_N_no1', s_N_no2 - s_N_no1



