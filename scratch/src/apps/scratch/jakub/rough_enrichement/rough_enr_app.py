'''
Created on Jun 19, 2009

@author: jakub
'''
if __name__ == '__main__':
    def example_1d():
        from ibvpy.api import FEDomain, FERefinementGrid, FEGrid, TStepper as TS,\
            BCDof, RTraceDomainListField
        from ibvpy.core.tloop import TLoop, TLine
        from ibvpy.mesh.xfe_subdomain import XFESubDomain
        from mats1D_elastic_tf import MATS1DElasticTF
        from ibvpy.mats.mats1D.mats1D_elastic.mats1D_elastic import MATS1DElastic
        from ibvpy.fets.fets1D.fets1D2l import FETS1D2L
        from ibvpy.fets.fets1D.fets1D2l3u import FETS1D2L3U
        from fets_crack_rough import FETSCrackRough
        fets_eval = FETS1D2L(mats_eval = MATS1DElastic(E= 2.0))
        xfets_eval = FETSCrackRough( parent_fets = fets_eval, 
                                     mats_eval = MATS1DElasticTF(E_m= 1.,E_f = 1., G = 1.))
    
        # Discretization
        
        fe_domain = FEDomain()
        fe_level1 = FERefinementGrid( domain = fe_domain, fets_eval = fets_eval )
        fe_grid1 = FEGrid( coord_max = (3.,0.,0.), 
                           shape   = (3,),
                           fets_eval = fets_eval,
                           level = fe_level1 )
        
        enr = True
        if enr:
            fe_xdomain = XFESubDomain( domain = fe_domain,
                                       fets_eval = xfets_eval,
                                       #fe_grid_idx_slice = fe_grid1[1,0],
                                       fe_grid_slice = fe_grid1['X  - 1.5'] )

        ts = TS( dof_resultants = True,
                 sdomain = fe_domain,
                 bcond_list =  [BCDof(var='u', value = 1., dims = [0],
                                      dof = 3 ),
                                BCDof(var='u', value = 0., dims = [0],
                                      dof = 0 ),
                                 BCDof(var='u', value = 0., dims = [0],
                                      dof = 4 ),
                                BCDof(var='u', value = 0., dims = [0],
                                      dof = 5 ),       
                                BCDof(var='u', value = 0., dims = [0],
                                      dof = 6 ),
                                BCDof(var='u', value = 0., dims = [0],
                                      dof = 7 ),   
#                                BCDof(var='u', value = 0., dims = [0],
#                                      link_dofs = [8,9,13],
#                                      link_coeffs = [-1,-1,-1],
#                                      dof = 12),
                                           ],
                 rtrace_list =  [ 
#                                 RTraceGraph(name = 'Fi,right over u_right (iteration)' ,
#                                       var_y = 'F_int', idx_y = 0,
#                                       var_x = 'U_k', idx_x = 1),
#                            RTraceDomainListField(name = 'Stress' ,
#                                 var = 'sig_app', idx = 0, warp = True ),
                             RTraceDomainListField(name = 'Displ m' ,
                                            var = 'u_rm', idx = 0,
                                            warp = True,
                                            warp_var = 'u_rm'),
                            RTraceDomainListField(name = 'Displ f' ,
                                            var = 'u_rf', idx = 0,
                                            warp = True,
                                            warp_var = 'u_rf'),              
                            
#                                     RTraceDomainField(name = 'N0' ,
#                                                  var = 'N_mtx', idx = 0,
#                                                  record_on = 'update')
                        ]             
                    )
#        
#        # Add the time-loop control
        tloop = TLoop( tstepper = ts,
                       tolerance = 1e-4, 
                       #KMAX = 2,
                       #debug = True, RESETMAX = 2,
                       tline  = TLine( min = 0.0,  step = 1, max = 1.0 ))
        
        #print "elements ",fe_xdomain.elements[0]
        if enr:
            fe_xdomain.deactivate_sliced_elems()
            print 'parent elems ',fe_xdomain.fe_grid_slice.elems
            print 'parent dofs ',fe_xdomain.fe_grid_slice.dofs
            print "dofmap ",fe_xdomain.elem_dof_map
            print "ls_values ", fe_xdomain.dots.dof_node_ls_values
            print 'intersection points ',fe_xdomain.fe_grid_slice.r_i#
            print "triangles ", fe_xdomain.dots.int_division
            print 'ip_coords', fe_xdomain.dots.ip_coords
            print 'ip_weigths', fe_xdomain.dots.ip_weights
            print 'ip_offset ',fe_xdomain.dots.ip_offset
            print 'ip_X_coords', fe_xdomain.dots.ip_X
            print 'ip_ls', fe_xdomain.dots.ip_ls_values
            print 'vtk_X ', fe_xdomain.dots.vtk_X
            print 'vtk triangles ', fe_xdomain.dots.rt_triangles
            print "vtk data ", fe_xdomain.dots.get_vtk_cell_data('blabla',0,0)
            print 'vtk_ls', fe_xdomain.dots.vtk_ls_values
            print 'J_det ',fe_xdomain.dots.J_det_grid
        
        print tloop.eval()
#        #ts.setup()
        from ibvpy.plugins.ibvpy_app import IBVPyApp
        ibvpy_app = IBVPyApp( ibv_resource = ts )
        ibvpy_app.main()
#            
#    #    print ts.F_int
#    #    print ts.rtrace_list[0].trace.ydata
    example_1d()