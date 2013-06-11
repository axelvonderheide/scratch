'''
Created on Mar 29, 2011

@author: jakub
'''
if __name__ == '__main__':
    def example_2d():
        from ibvpy.api import FEDomain, FERefinementGrid, FEGrid, TStepper as TS, \
            BCDofGroup, RTraceDomainListField, BCDof
        from ibvpy.core.tloop import TLoop, TLine
        from ibvpy.mesh.xfe_subdomain import XFESubDomain
        from ibvpy.mats.mats2D.mats2D_elastic.mats2D_elastic import MATS2DElastic
        from ibvpy.fets.fets2D.fets2D4q import FETS2D4Q
        from fets_crack_tf import FETSCrackTF
        from ibvpy.mats.mats2D5.mats2D5_bond.mats2D_bond import MATS2D5Bond
        from ibvpy.mats.mats2D5.mats2D5_bond.mats2D5_plastic_bond import MATS2D5PlasticBond

        from ibvpy.fets.fets2D.fets2Dtf import FETS2DTF
        from ibvpy.fets.fets2D.fets2D4q import FETS2D4Q
        from ibvpy.fets.fets2D.fets2D4q8u import FETS2D4Q8U
        from ibvpy.fets.fets2D.fets2D4q9u import FETS2D4Q9U
        from ibvpy.fets.fets2D.fets2D4q12u import FETS2D4Q12U
        from ibvpy.fets.fets2D.fets2D4q16u import FETS2D4Q16U
#        from fets2D4qtf import FETS2D4QTF
#        from fets2D4q8utf import FETS2D4Q8UTF
#        from fets2D4q9utf import FETS2D4Q9UTF
        #fets_eval = FETS2D4Q(mats_eval = MATS2DElastic(E= 1.,nu=0.))
        fets_eval = FETS2DTF( parent_fets = FETS2D4Q16U(),
                              mats_eval = MATS2D5PlasticBond( 
                                                      E_m = 340., nu_m = 0.2,
                                                       E_f = 6., nu_f = 0.,
                                                       G = 1.25e5, sigma_y = 100000.,
                                                       stress_state = 'plane_strain' ) )
        xfets_eval = FETSCrackTF( parent_fets = fets_eval,
                                  debug_on = True,
                                  int_order = 5 )

        # Discretization

        fe_domain = FEDomain()
        fe_level1 = FERefinementGrid( domain = fe_domain,
                                      fets_eval = fets_eval )
        fe_grid1 = FEGrid( coord_min = ( 1.5, -1.5, 0. ),
                          coord_max = ( 4.5, 1.5, 0. ),
                           shape = ( 3, 3 ),
                           fets_eval = fets_eval,
                           level = fe_level1 )

        bcond_lock = []
        enr = True
        if enr:
            fe_xdomain = XFESubDomain( domain = fe_domain,
                                       rt_quad = False,
                                       fets_eval = xfets_eval,
                                       #fe_grid_idx_slice = fe_grid1[1,0],
                                       #fe_grid_slice = fe_grid1['X**2 + Y**2  - 9.1 '] 
                                       fe_grid_slice = fe_grid1['X - 3. ']
                                       )
            for i in range( 400, 480 ):
                bcond_lock.append( BCDof( var = 'u', dof = i, value = 0. ) )


        ts = TS( dof_resultants = True,
                 sdomain = fe_domain,
                 bcond_list = [
                                BCDofGroup( var = 'u', value = .2, dims = [0, 2 ],
                                          get_dof_method = fe_grid1.get_right_dofs ),
                                BCDofGroup( var = 'u', value = 0., dims = [1, 3],
                                          get_dof_method = fe_grid1.get_right_dofs ),
                                BCDofGroup( var = 'u', value = 0., dims = [0, 1, 2, 3],
                                           get_dof_method = fe_grid1.get_left_dofs ),
                                           ] + bcond_lock,
                 rtrace_list = [
#                                 RTraceGraph(name = 'Fi,right over u_right (iteration)' ,
#                                       var_y = 'F_int', idx_y = 0,
#                                       var_x = 'U_k', idx_x = 1),
#                            RTraceDomainListField(name = 'Stress' ,
#                                 var = 'sig_app', idx = 0, warp = True ),
                             RTraceDomainListField( name = 'Displ matrix' ,
                                            var = 'u_m',
                                            warp = True,
                                            warp_var = 'u_m' ),
                            RTraceDomainListField( name = 'Displ reinf' ,
                                            var = 'u_f',
                                            warp = True,
                                            warp_var = 'u_f' ),
                            RTraceDomainListField( name = 'Strain matrix' ,
                                            var = 'eps_m',
                                            warp = True,
                                            warp_var = 'u_m' ),
                            RTraceDomainListField( name = 'Strain reinf' ,
                                            var = 'eps_f',
                                            warp = True,
                                            warp_var = 'u_f' ),
                            RTraceDomainListField( name = 'slip' ,
                                            var = 'eps_b',
                                            #warp = True,
                                            warp_var = 'u_f' ),

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
                       tline = TLine( min = 0.0, step = 1, max = 1.0 ) )

        if enr:
            #print "elements ",fe_xdomain.elements[0]
            fe_xdomain.deactivate_sliced_elems()
            #fe_xdomain2.deactivate_sliced_elems()
            print 'parent elems ', fe_xdomain.fe_grid_slice.elems
            print 'parent dofs ', fe_xdomain.fe_grid_slice.dofs
            print "dofmap ", fe_xdomain.elem_dof_map
            print "ls_values ", fe_xdomain.dots.dof_node_ls_values
            print 'intersection points ', fe_xdomain.fe_grid_slice.r_i
            print "triangles ", fe_xdomain.dots.rt_triangles
            print "vtk points ", fe_xdomain.dots.vtk_X
            print "vtk data ", fe_xdomain.dots.get_vtk_cell_data( 'nodes', 0, 0 )
            print 'ip_triangles', fe_xdomain.dots.int_division
            print 'ip_coords', fe_xdomain.dots.ip_coords
            print 'ip_weigths', fe_xdomain.dots.ip_weights
            print 'ip_offset', fe_xdomain.dots.ip_offset
            print 'ip_X_coords', fe_xdomain.dots.ip_X
            print 'ip_ls', fe_xdomain.dots.ip_ls_values
            print 'vtk_ls', fe_xdomain.dots.vtk_ls_values
            print 'J_det ', fe_xdomain.dots.J_det_grid

        tloop.eval()
#        #ts.setup()
        from ibvpy.plugins.ibvpy_app import IBVPyApp
        ibvpy_app = IBVPyApp( ibv_resource = ts )
        ibvpy_app.main()
#            
#    #    print ts.F_int
#    #    print ts.rtrace_list[0].trace.ydata

    example_2d()
