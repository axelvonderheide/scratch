'''
Created on Jul 2, 2009

@author: jakub
'''
if __name__ == '__main__':        
    def example_2d():
        from ibvpy.api import FEDomain, FERefinementGrid, FEGrid, TStepper as TS,\
            BCDofGroup, RTraceDomainListField
        from ibvpy.core.tloop import TLoop, TLine
        from ibvpy.mesh.xfe_subdomain import XFESubDomain
        from ibvpy.mats.mats2D.mats2D_elastic.mats2D_elastic import MATS2DElastic
        from ibvpy.fets.fets2D.fets2D4q import FETS2D4Q
        from fets_crack_tf import FETSCrackTF
        from mats2D_elastic_tf import MATS2DElasticTF
        from fets2D4qtf import FETS2D4QTF
        from fets2D4q8utf import FETS2D4Q8UTF
        from fets2D4q9utf import FETS2D4Q9UTF
        #fets_eval = FETS2D4Q(mats_eval = MATS2DElastic(E= 1.,nu=0.))
        fets_eval = FETS2D4QTF(mats_eval = MATS2DElasticTF(E_m = 30, nu_m = 0.,
                                                       E_f = 30, nu_f = 0.,
                                                       G = 8.)) 
        xfets_eval = FETSCrackTF( parent_fets = fets_eval, int_order = 5 )
    
        # Discretization
        
        fe_domain = FEDomain()
        fe_level1 = FERefinementGrid( domain = fe_domain, fets_eval = fets_eval )
        fe_grid1 = FEGrid( coord_max = (3.,2.,0.), 
                           shape   = (6,4),
                           fets_eval = fets_eval,
                           level = fe_level1 )
#        fe_grid1.deactivate( (1,0) )
#        fe_grid1.deactivate( (1,1) )
        
        enr = True
        if enr:
            fe_xdomain = XFESubDomain( domain = fe_domain,
                                       fets_eval = xfets_eval,
                                       #fe_grid_idx_slice = fe_grid1[1,0],
                                       fe_grid_slice = fe_grid1['X - 0.25 * Y - 1.251'] )

        ts = TS( dof_resultants = True,
                 sdomain = fe_domain,
                 bcond_list =  [BCDofGroup(var='u', value = .2, dims = [0],
                                          get_dof_method = fe_grid1.get_right_dofs ),
                                BCDofGroup(var='u', value = 0., dims = [1],
                                          get_dof_method = fe_grid1.get_bottom_right_dofs ),
                                BCDofGroup(var='u', value = 0., dims = [0],
                                           get_dof_method = fe_grid1.get_left_dofs ),
                                BCDofGroup(var='u', value = 0., dims = [1],
                                           get_dof_method = fe_grid1.get_bottom_left_dofs ),
                                           ],
                 rtrace_list =  [ 
#                                 RTraceGraph(name = 'Fi,right over u_right (iteration)' ,
#                                       var_y = 'F_int', idx_y = 0,
#                                       var_x = 'U_k', idx_x = 1),
#                            RTraceDomainListField(name = 'Stress' ,
#                                 var = 'sig_app', idx = 0, warp = True ),
                             RTraceDomainListField(name = 'Displ matrix' ,
                                            var = 'u_m',
                                            warp = True,
                                            warp_var = 'u_m'),
                            RTraceDomainListField(name = 'Displ reinf' ,
                                            var = 'u_f',
                                            warp = True,
                                            warp_var = 'u_f'),
                            RTraceDomainListField(name = 'Strain matrix' ,
                                            var = 'eps_m',
                                            warp = True,
                                            warp_var = 'u_f'),  
                        RTraceDomainListField(name = 'Strain reinf' ,
                                            var = 'eps_f',
                                            warp = True,
                                            warp_var = 'u_f'),    
                                
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
        if enr:
            #print "elements ",fe_xdomain.elements[0]
            fe_xdomain.deactivate_sliced_elems()
            print 'parent elems ',fe_xdomain.fe_grid_slice.elems
            print 'parent dofs ',fe_xdomain.fe_grid_slice.dofs
            print "dofmap ",fe_xdomain.elem_dof_map
            print "ls_values ", fe_xdomain.dots.dof_node_ls_values
            print 'intersection points ',fe_xdomain.fe_grid_slice.r_i
            print "triangles ", fe_xdomain.dots.rt_triangles
            print "vtk points ", fe_xdomain.dots.vtk_X
            print "vtk data ", fe_xdomain.dots.get_vtk_cell_data('blabla',0,0)
            print 'ip_triangles', fe_xdomain.dots.int_division
            print 'ip_coords', fe_xdomain.dots.ip_coords
            print 'ip_weigths', fe_xdomain.dots.ip_weights
            print 'ip_offset', fe_xdomain.dots.ip_offset
            print 'ip_X_coords', fe_xdomain.dots.ip_X
            print 'ip_ls', fe_xdomain.dots.ip_ls_values
            print 'vtk_ls', fe_xdomain.dots.vtk_ls_values
            print 'J_det ',fe_xdomain.dots.J_det_grid
        
        tloop.eval()
#        #ts.setup()
        from ibvpy.plugins.ibvpy_app import IBVPyApp
        ibvpy_app = IBVPyApp( ibv_resource = ts )
        ibvpy_app.main()
#            
#    #    print ts.F_int
#    #    print ts.rtrace_list[0].trace.ydata

    example_2d()