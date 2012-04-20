'''
Created on Sep 23, 2009

@author: jakub
'''
if __name__ == '__main__':
        from ibvpy.api import FEDomain, FERefinementGrid, FEGrid, TStepper as TS, \
            BCDofGroup, BCDof, RTraceDomainListField
        from ibvpy.core.tloop import TLoop, TLine
        from ibvpy.mesh.xfe_subdomain import XFESubDomain
        from ibvpy.mats.mats2D.mats2D_elastic.mats2D_elastic import MATS2DElastic
        from ibvpy.fets.fets2D.fets2D4q import FETS2D4Q
        from ibvpy.fets.fets2D.fets2D4q8u import FETS2D4Q8U
        from ibvpy.fets.fets2D.fets2D4q9u import FETS2D4Q9U
        from ibvpy.fets.fets2D.fets2D9q import FETS2D9Q
        from ibvpy.fets.fets_ls.fets_crack import FETSCrack

        fets_eval = FETS2D4Q( mats_eval = MATS2DElastic( E = 1., nu = 0. ) )
        xfets_eval = FETSCrack( parent_fets = fets_eval, int_order = 5 )

        # Discretization

        fe_domain = FEDomain()
        fe_level1 = FERefinementGrid( domain = fe_domain, fets_eval = fets_eval )
        fe_grid1 = FEGrid( coord_max = ( 1., 1., 0. ),
                           shape = ( 1, 3 ),
                           fets_eval = fets_eval,
                           level = fe_level1 )
#        fe_grid1.deactivate( (1,0) )
#        fe_grid1.deactivate( (1,1) )

        fe_xdomain = XFESubDomain( domain = fe_domain,
                                   fets_eval = xfets_eval,
                                   boundary = lambda X, Y: Y < 0.5,
                                   #fe_grid_idx_slice = fe_grid1[1,0],
                                   slice = fe_grid1['X  -  0.5 '] )

        ts = TS( dof_resultants = True,
                 sdomain = fe_domain,
                 bcond_list = [BCDofGroup( var = 'u', value = 0., dims = [0, 1],
                                          get_dof_method = fe_grid1.get_right_dofs ),
                                BCDofGroup( var = 'u', value = 0., dims = [1],
                                          get_dof_method = fe_grid1.get_left_dofs ),
                                BCDofGroup( var = 'u', value = -.2, dims = [0],
                                           get_dof_method = fe_grid1.get_left_dofs ),
                                           ],
                 rtrace_list = [
#                                 RTraceGraph(name = 'Fi,right over u_right (iteration)' ,
#                                       var_y = 'F_int', idx_y = 0,
#                                       var_x = 'U_k', idx_x = 1),
#                            RTraceDomainListField(name = 'Stress' ,
#                                 var = 'sig_app', idx = 0, warp = True ),
                             RTraceDomainListField( name = 'Displacement' ,
                                            var = 'u', idx = 0,
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
#        print 'parent elems ', fe_xdomain.fe_grid_slice.elems
#        print 'parent dofs ', fe_xdomain.fe_grid_slice.dofs
#        print "dofmap ", fe_xdomain.elem_dof_map
#        print "ls_values ", fe_xdomain.dots.dof_node_ls_values
#        print 'intersection points ', fe_xdomain.fe_grid_slice.r_i
#        print "triangles ", fe_xdomain.dots.rt_triangles
#        print "vtk points ", fe_xdomain.dots.vtk_X
#        #print "vtk data ", fe_xdomain.dots.get_vtk_cell_data( 'blabla', 0, 0 )
#        print 'ip_triangles', fe_xdomain.dots.int_division
#        print 'ip_coords', fe_xdomain.dots.ip_coords
#        print 'ip_weigths', fe_xdomain.dots.ip_weights
#        print 'ip_offset', fe_xdomain.dots.ip_offset
#        print 'ip_X_coords', fe_xdomain.dots.ip_X
#        print 'ip_ls', fe_xdomain.dots.ip_ls_values
#        print 'vtk_ls', fe_xdomain.dots.vtk_ls_values
#        print 'J_det ', fe_xdomain.dots.J_det_grid

        print tloop.eval()
#        #ts.setup()
        from ibvpy.plugins.ibvpy_app import IBVPyApp
        ibvpy_app = IBVPyApp( ibv_resource = ts )
        ibvpy_app.main()
