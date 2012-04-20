'''
Created on Oct 22, 2009

@author: jakub
'''
if __name__ == '__main__':
        from numpy import cos, sin, arange, vstack, array, pi
        from ibvpy.api import FEDomain, FERefinementGrid, FEGrid, TStepper as TS, \
            BCDofGroup, BCSlice
        from apps.scratch.jakub.mlab.mlab_trace import RTraceDomainListField
        from ibvpy.core.tloop import TLoop, TLine
        from ibvpy.mesh.xfe_subdomain import XFESubDomain
        from ibvpy.mats.mats2D.mats2D_elastic.mats2D_elastic import MATS2DElastic
        from ibvpy.fets.fets2D.fets2D4q import FETS2D4Q
        from ibvpy.fets.fets2D.fets2D4q8u import FETS2D4Q8U
        from ibvpy.fets.fets2D.fets2D4q9u import FETS2D4Q9U

        from ibvpy.fets.fets_ls.fets_crack import FETSCrack#xfem jump
        from ibvpy.fets.fets_ls.fets_ls_eval import FETSLSEval#selective ingration

        fets_eval = FETS2D4Q( mats_eval = MATS2DElastic( E = 200000000., nu = 0. ) )
        xfets_eval = FETSLSEval( parent_fets = fets_eval, int_order = 3 )#int order =1,3,5
        xfets_eval2 = FETSLSEval( parent_fets = fets_eval, int_order = 3 )

        # Discretization

        fe_domain = FEDomain()
        fe_level1 = FERefinementGrid( domain = fe_domain, fets_eval = fets_eval )
        y_fines = 22
        height = 2.0
        fe_grid1 = FEGrid( coord_max = ( 2., 2., 0. ),
                           n_elems = ( 22, y_fines ),
                           fets_eval = fets_eval,
                           level = fe_level1 )

        fe_xdomain = XFESubDomain( domain = fe_domain,
                                   domain_tag = 'pos',
                                   rt_tol = 0.000,
                                   fets_eval = xfets_eval,
                                   fe_grid_slice = fe_grid1['((X-0.7)**2+((Y-0.6))**2)**0.5-.4' ] )

        sektors = 6.0
        a = arange( 0, sektors + 1 )
        b = sin( a / sektors * 2 * pi ) * 0.4 + 1.4
        c = cos( a / sektors * 2 * pi ) * 0.2 + 1.5
        d = array( vstack( ( b, c ) ).T )
        print "d", d

        fe_xdomain2 = XFESubDomain( domain = fe_domain,
                                    domain_tag = 'pos',
                                    rt_tol = 0.0001,
                                    fets_eval = xfets_eval2,
                                    fe_grid_slice = fe_grid1[d] )

        ts = TS( dof_resultants = True,
                 sdomain = fe_domain,
                 bcond_list = [BCDofGroup( var = 'u', value = 0., dims = [1], #var u or f,dims [0,1]
                                          get_dof_method = fe_grid1.get_right_dofs ),
                                BCSlice( var = 'f', value = 100000. / y_fines * height / 2, dims = [0],
                                        slice = fe_grid1[-1, :, -1, :] ),
#                                 BCSlice(var='f', value=0.1,dims=[0]
#                                        slice=fe_grid1[-1,0,-1,0]),        
#                                 BCSlice(var='f', value=0.1,dims=[0],
#                                        slice=fe_grid1[-1,-1,-1,-1]),
#                                 BCSlice(var='f', value=0.1,dims=[0],
#                                        slice=fe_grid1[-1,0,-1,1]),        
#                                 BCSlice(var='f', value=0.1,dims=[0],
#                                        slice=fe_grid1[-1,-1,-1,0]),
                                BCDofGroup( var = 'u', value = 0., dims = [1],
                                          get_dof_method = fe_grid1.get_bottom_left_dofs ),
                                BCDofGroup( var = 'u', value = 0., dims = [0],
                                           get_dof_method = fe_grid1.get_left_dofs )
                                #BCSlice(var='u', value = 0., dims = [0,1],
                                 #       slice = fe_grid1[1,0,:,0] )#slice = elem_x, elemy, nodex,nodey       
                                           ],
                 rtrace_list = [
#                                 RTraceGraph(name = 'Fi,right over u_right (iteration)' ,
#                                       var_y = 'F_int', idx_y = 0,
#                                       var_x = 'U_k', idx_x = 1),
                            RTraceDomainListField( name = 'Stress' , #name = arbitrary string
                                 var = 'sig_app', warp = True ), #var = sig_app, eps, u
                             RTraceDomainListField( name = 'Displacement' ,
                                            var = 'u',
                                            warp = True ),
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
                       tline = TLine( min = 0.0, step = 1., max = 1.0 ) )

        #print "elements ",fe_xdomain.elements[0]
        fe_xdomain.deactivate_sliced_elems()
#        print 'ls string ',fe_xdomain.fe_grid_slice.ls_function_eval
#        print 'parent elems ',fe_xdomain.fe_grid_slice.elems
#        print 'parent dofs ',fe_xdomain.fe_grid_slice.dofs
#        print "dofmap ",fe_xdomain.elem_dof_map
#        print "ls_values ", fe_xdomain.dots.dof_node_ls_values
#        print 'intersection points ',fe_xdomain.fe_grid_slice.r_i
#        print "triangles ", fe_xdomain.dots.rt_triangles
#        print "vtk points ", fe_xdomain.dots.vtk_X
#        print "vtk data ", fe_xdomain.dots.get_vtk_cell_data('blabla',0,0)
#        print 'ip_triangles', fe_xdomain.dots.int_division
#        print 'ip_coords', fe_xdomain.dots.ip_coords
#        print 'ip_weigths', fe_xdomain.dots.ip_weights
#        print 'ip_offset', fe_xdomain.dots.ip_offset
#        print 'ip_X_coords', fe_xdomain.dots.ip_X
#        print 'ip_ls', fe_xdomain.dots.ip_ls_values
#        print 'vtk_ls', fe_xdomain.dots.vtk_ls_values
#        print 'J_det ',fe_xdomain.dots.J_det_grid

        fe_xdomain2.deactivate_sliced_elems()
#        print 'ls string ',fe_xdomain2.fe_grid_slice.ls_function_eval
#        print '2parent elems ',fe_xdomain2.fe_grid_slice.elems
#        print '2parent dofs ',fe_xdomain2.fe_grid_slice.dofs
#        print "2dofmap ",fe_xdomain2.elem_dof_map
#        print "2ls_values ", fe_xdomain2.dots.dof_node_ls_values
#        print '2intersection points ',fe_xdomain2.fe_grid_slice.r_i
#        print "2triangles ", fe_xdomain2.dots.rt_triangles
#        print "2vtk points ", fe_xdomain2.dots.vtk_X
#        print "2vtk data ", fe_xdomain2.dots.get_vtk_cell_data('blabla',0,0)
#        print '2ip_triangles', fe_xdomain2.dots.int_division
#        print '2ip_coords', fe_xdomain2.dots.ip_coords
#        print '2ip_weigths', fe_xdomain2.dots.ip_weights
#        print '2ip_offset', fe_xdomain2.dots.ip_offset
#        print '2ip_X_coords', fe_xdomain2.dots.ip_X
#        print '2ip_ls', fe_xdomain2.dots.ip_ls_values
#        print '2vtk_ls', fe_xdomain2.dots.vtk_ls_values
#        print '2J_det ',fe_xdomain2.dots.J_det_grid

        tloop.eval()
#        #ts.setup()
        from ibvpy.plugins.ibvpy_app import IBVPyApp
        ibvpy_app = IBVPyApp( ibv_resource = ts )
        ibvpy_app.main()
