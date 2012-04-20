'''
Created on Nov 12, 2009

@author: jakub
'''
if __name__ == '__main__':        
    def example_2d():
        from ibvpy.api import FEDomain, FERefinementGrid, FEGrid, TStepper as TS,\
            BCDofGroup, RTraceDomainListField, BCSlice
        from ibvpy.core.tloop import TLoop, TLine
        from ibvpy.mesh.xfe_subdomain import XFESubDomain
        from ibvpy.mats.mats2D.mats2D_elastic.mats2D_elastic import MATS2DElastic
        from ibvpy.fets.fets2D.fets2D4q import FETS2D4Q
        from fets_crack_tf import FETSCrackTF
        from ibvpy.mats.mats2D5.mats2D5_bond.mats2D_bond import MATS2D5Bond
        from ibvpy.fets.fets2D.fets2Dtf import FETS2DTF
        from ibvpy.fets.fets2D.fets2D4q import FETS2D4Q
        from ibvpy.fets.fets2D.fets2D4q8u import FETS2D4Q8U
        from ibvpy.fets.fets2D.fets2D4q9u import FETS2D4Q9U
        from ibvpy.fets.fets2D.fets2D4q12u import FETS2D4Q12U
        from ibvpy.fets.fets2D.fets2D4q16u import FETS2D4Q16U
        from ibvpy.fets.fets_eval import RTraceEvalElemFieldVar
        from numpy import dot, array, linalg, sum, deg2rad
        from math import pi, sin, cos
#        from fets2D4qtf import FETS2D4QTF
#        from fets2D4q8utf import FETS2D4Q8UTF
#        from fets2D4q9utf import FETS2D4Q9UTF
        #fets_eval = FETS2D4Q(mats_eval = MATS2DElastic(E= 1.,nu=0.))
        fets_eval = FETS2DTF(parent_fets = FETS2D4Q12U(),
                              mats_eval = MATS2D5Bond(E_m = 2., nu_m = 0.,
                                                       E_f = 1., nu_f = 0.,
                                                       G = 10., 
                                                       stress_state = 'plane_stress')) 
        xfets_eval = FETSCrackTF( parent_fets = fets_eval, int_order = 5 )
    
        alpha = 90.# angle of the crack measured anticlockwise from x axis
        alpha_r = deg2rad(alpha)
        crack_normal = array([sin(alpha_r), -cos(alpha_r)])
        print 'crack_normal ',crack_normal
    
        def xget_eps_mn(sctx, u):
            e_id = sctx.e_id
            p_id = sctx.p_id
            X_mtx = sctx.X
            r_pnt = sctx.loc
            B_mtx = xfets_eval.get_B_mtx(r_pnt, X_mtx,
                                   sctx.dots.dof_node_ls_values[e_id],
                                   sctx.dots.vtk_ls_values[e_id][p_id] )
            eps = dot( B_mtx, u )
            p_eps, p_vct = linalg.eigh(array([[eps[0],eps[2]],[eps[2],eps[1]]]))
            return array([-sum(dot((p_eps[:,None] * p_vct),crack_normal))])
    
        def xget_eps_fn(sctx, u):
            e_id = sctx.e_id
            p_id = sctx.p_id
            X_mtx = sctx.X
            r_pnt = sctx.loc
            B_mtx = xfets_eval.get_B_mtx(r_pnt, X_mtx,
                                   sctx.dots.dof_node_ls_values[e_id],
                                   sctx.dots.vtk_ls_values[e_id][p_id] )
            eps = dot( B_mtx, u )
            p_eps, p_vct = linalg.eigh(array([[eps[3],eps[5]],[eps[5],eps[4]]]))
            return array([-sum(dot((p_eps[:,None] * p_vct),crack_normal))])
    
        xfets_eval.rte_dict.update( {'eps_mn' : RTraceEvalElemFieldVar( eval = xget_eps_mn ), 
                                     'eps_fn' : RTraceEvalElemFieldVar( eval = xget_eps_fn )} )
            
        def get_eps_mn( sctx, u):
            X_mtx = sctx.X
            r_pnt = sctx.loc
            B_mtx = fets_eval.get_B_mtx(r_pnt, X_mtx)
            eps = dot( B_mtx, u )
            p_eps, p_vct = linalg.eigh(array([[eps[0],eps[2]],[eps[2],eps[1]]]))
            return array([-sum(dot((p_eps[:,None] * p_vct),crack_normal))])
        
        def get_eps_fn( sctx, u):
            X_mtx = sctx.X
            r_pnt = sctx.loc
            B_mtx = fets_eval.get_B_mtx(r_pnt, X_mtx)
            eps = dot( B_mtx, u )
            p_eps, p_vct = linalg.eigh(array([[eps[3],eps[5]],[eps[5],eps[4]]]))
            return array([-sum(dot((p_eps[:,None] * p_vct),crack_normal))])
        
        fets_eval.rte_dict.update( {'eps_mn' : RTraceEvalElemFieldVar( eval = get_eps_mn ), 
                                    'eps_fn' : RTraceEvalElemFieldVar( eval = get_eps_fn )} )
    
        # Discretization
        
        fineness_x = 15
        fineness_y = 10
        
        fe_domain = FEDomain()
        fe_level1 = FERefinementGrid( domain = fe_domain, fets_eval = fets_eval )
        fe_grid1 = FEGrid( coord_max = (3.1,2.1,0.), 
                           shape   = (fineness_x,fineness_y),
                           fets_eval = fets_eval,
                           level = fe_level1 )
        
        enr = True
        if enr:
            fe_xdomain = XFESubDomain( domain = fe_domain,
                                       fets_eval = xfets_eval,
                                       fe_grid_slice = fe_grid1['X - 1.55'] )#90.deg
                                       #fe_grid_slice = fe_grid1['X - 1.7320508075688767 *Y +0.26865334794732054'] )#30.deg
                                       #fe_grid_slice = fe_grid1['Y - 1.05'] )#0.deg
                                       #fe_grid_slice = fe_grid1['X-Y - .55'] )#45.deg
                                       #fe_grid_slice = fe_grid1['X-0.57735026918962573*Y - 0.94378221735089296'] )#60.deg
        ts = TS( dof_resultants = True,
                 sdomain = fe_domain,
                 bcond_list =  [
#                                BCDofGroup(var='u', value = .2, dims = [0],
#                                          get_dof_method = fe_grid1.get_right_dofs ),
                                #force bc valid for bilinear elems
#                                BCSlice(var='f', value = .1/fineness_y, dims = [2],
#                                          slice =  fe_grid1[-1,:,-1,:]),
                                ###################################################
                                #force bc valid for bicubic elems
                                BCSlice(var='f', value = .025/2./fineness_y, dims = [2],
                                          slice =  fe_grid1[-1,:,-1,0]),
                                BCSlice(var='f', value = .025/2./fineness_y, dims = [2],
                                          slice =  fe_grid1[-1,:,-1,-1]),
                                BCSlice(var='f', value = .075/2./fineness_y, dims = [2],
                                          slice =  fe_grid1[-1,:,-1,1:-1]),
                                #####################################################
                                BCDofGroup(var='u', value = 0., dims = [2],
                                           get_dof_method = fe_grid1.get_left_dofs ),
                                BCDofGroup(var='u', value = 0., dims = [3],
                                           get_dof_method = fe_grid1.get_bottom_left_dofs ),
                                BCDofGroup(var='u', value = 0., dims = [3],
                                           get_dof_method = fe_grid1.get_bottom_right_dofs ),
                                           ],
                 rtrace_list =  [ 
#                                 RTraceGraph(name = 'Fi,right over u_right (iteration)' ,
#                                       var_y = 'F_int', idx_y = 0,
#                                       var_x = 'U_k', idx_x = 1),
#                            RTraceDomainListField(name = 'Stress' ,
#                                 var = 'sig_app', idx = 0, warp = True ),
#                             RTraceDomainListField(name = 'Displ matrix' ,
#                                            var = 'u_m',
#                                            warp = True,
#                                            warp_var = 'u_m'),
#                            RTraceDomainListField(name = 'Displ reinf' ,
#                                            var = 'u_f',
#                                            warp = True,
#                                            warp_var = 'u_f'),
#                            RTraceDomainListField(name = 'Strain matrix' ,
#                                            var = 'eps_mn',
#                                            warp = True,
#                                            warp_var = 'u_m'),  
                            RTraceDomainListField(name = 'Strain reinf' ,
                                            var = 'eps_fn')#,
#                                            warp = True,
#                                            warp_var = 'u_f'),    
#                            RTraceDomainListField(name = 'Strain reinf mtrl' ,
#                                            var = 'eps_app_f',
#                                            warp = True,
#                                            warp_var = 'u_f'), 
#                                     RTraceDomainField(name = 'N0' ,
#                                                  var = 'N_mtx', idx = 0,
#                                                  record_on = 'update')
                        ]             
                    )
#        
#        # Add the time-loop control
        tloop = TLoop( tstepper = ts,
#                       tolerance = 1e-4, KMAX = 4,
#                       debug = True, RESETMAX = 2,
                       tline  = TLine( min = 0.0,  step = 1, max = 1.0 ))
        fe_xdomain.deactivate_sliced_elems()
        
        tloop.eval()
#        #ts.setup()
        from ibvpy.plugins.ibvpy_app import IBVPyApp
        ibvpy_app = IBVPyApp( ibv_resource = ts )
        ibvpy_app.main()
#            
#    #    print ts.F_int
#    #    print ts.rtrace_list[0].trace.ydata

    example_2d()