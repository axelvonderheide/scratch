
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
                                                   var_x = 'e_equiv_arr', idx_x = 0,
                                                   var_y = 's_equiv_arr', idx_y = 0,
                                                   record_on = 'update' ),

                                       # microplane fracture energy:
                                       MATS2DMicroplaneDamageTraceGfmic(name = 'time - G_f_micro',
                                                                     var_x = 'e_equiv_arr', idx_x = 0,
                                                                     var_y = 's_equiv_arr', idx_y = 0,
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
#                                       RTraceArraySnapshot(name = 'fracture energy contributions',
#                                                           var = 'fracture_energy_list',
#                                                           record_on = 'update' ),
#                                       RTraceArraySnapshot(name = 'microplane damage',
#                                                           var = 'microplane_damage',
#                                                           record_on = 'update' ),
                                     ]
                       )

# --- test if product-type symmetrization works properly in the elastic regime ---    
    from numpy import \
        array, ones, zeros, outer, inner, transpose, dot, frompyfunc, \
        fabs, linspace, vdot, identity, tensordot, \
        sin as nsin, meshgrid, float_, ix_, \
        vstack, hstack, sqrt as arr_sqrt, swapaxes, copy
    from scipy.linalg import eigh, inv
    
    mats2D_explore.mats2D_eval.polar_fn.n_mp = 7
    mats2D_explore.mats2D_eval.elastic_debug = False
    mats2D_explore.mats2D_eval.symmetrization = 'product-type'
    mats2D_explore.mats2D_eval.model_version = 'compliance'

    mats2D_explore.tloop.eval()
    mats2D_explore.tloop.setup()
    
    phi = zeros((2,2),dtype=float)
    phi[0,0] = 1.
    phi[0,1] = 1.e-15
    
    # without that line phi becomes unsymmetric!!!
#    phi[1,0] = phi[0,1]

    # an unsymmetric matrix leads to unreasonable results if eig() is used insted of eigh().
    # eigh() solves the symmetric eigenvalueproblem based on the lower triangle.
    # eig() solves the ordinary eigenvalue problem (needs about twice the time of eigh.    

#    phi[0,1] = 0.
    phi[1,1] = 1.
    print 'phi',phi
    
    psi = phi    
    
    # ---
    n_dim = 2
    phi_mtx = phi
 
    phi_eig_value, phi_eig_mtx = eigh( phi_mtx )
    print 'phi_eig_value, phi_eig_mtx', phi_eig_value, phi_eig_mtx
    
    phi_eig_value_real = array([ pe.real for pe in phi_eig_value] )                
    # transform phi_mtx to PDC:
    # (assure that besides the diagonal the entries are exactly zero)
    phi_pdc_mtx = zeros((n_dim,n_dim),dtype=float)
    for i in range(n_dim):
        phi_pdc_mtx[i,i] = phi_eig_value_real[i]
    print 'phi_pdc_mtx', phi_pdc_mtx
    # w_mtx = tensorial square root of the second order damage tensor:
    w_pdc_mtx = arr_sqrt( phi_pdc_mtx )
    print "w_pdc_mtx", w_pdc_mtx
    # transform the matrix w back to x-y-coordinates:
    w_mtx = dot(dot(phi_eig_mtx, w_pdc_mtx),transpose(phi_eig_mtx))
    print "w_mtx", w_mtx
    # ---
    
#    M4 = mats2D_explore.mats2D_eval._get_M_tns_sum_type(phi)
    M4 = mats2D_explore.mats2D_eval._get_M_tns_product_type(psi)
    print 'M4'
    print M4
    
    C4_e = mats2D_explore.mats2D_eval.C4_e
    print 'elastic tensor'
    print C4_e

    C4_mdm = tensordot( M4, tensordot( C4_e, M4, [[2,3],[0,1]] ), [[0,1],[0,1]] )
    
    print 'after multiplication'
    print C4_mdm

    # compare the elastic compliance matrix with the undamaged compliance matrix C4_mdm
    print 'difference'
    print C4_e - C4_mdm




