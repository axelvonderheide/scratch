'''
Created on Oct 22, 2009

@author: jakub
'''
def example_with_new_domain():
    from ibvpy.api import \
        TStepper as TS, RTraceGraph, RTraceDomainListField, \
        RTraceDomainListInteg, TLoop,FEDomain, FERefinementGrid,\
        TLine, BCDof, IBVPSolve as IS, DOTSEval
    from ibvpy.mats.mats2D.mats2D_conduction.mats2D_conduction import MATS2DConduction
    from ibvpy.mats.mats2D.mats2D_elastic.mats2D_elastic import MATS2DElastic
    from ibvpy.fets.fets2D.fets2D4q4t import FETS2D4Q4T
    from ibvpy.fets.fets2D.fets2D4q import FETS2D4Q
    from ibvpy.api import BCDofGroup
    from ibvpy.fets.fets_ls.fets_ls_eval import FETSLSEval
    from ibvpy.fets.fets_ls.fets_crack import FETSCrack
    from ibvpy.mesh.xfe_subdomain import XFESubDomain
    
    #fets_eval = FETS2D4Q(mats_eval = MATS2DElastic(E= 1.,nu=0.))
    fets_eval = FETS2D4Q4T(mats_eval = MATS2DConduction(k = 1.)) 
    xfets_eval = FETSLSEval( parent_fets = fets_eval, int_order = 3 )

    
    from ibvpy.mesh.fe_grid import FEGrid
    from ibvpy.mesh.fe_refinement_grid import FERefinementGrid
    from ibvpy.mesh.fe_domain import FEDomain
    from mathkit.mfn.mfn_line.mfn_line import MFnLineArray

    # Discretization
    fe_domain = FEDomain()
    fe_level1 = FERefinementGrid( domain = fe_domain, fets_eval = fets_eval )
    fe_grid = FEGrid( coord_max = (1.,1.,0.), 
                       shape   = (1,1),
                       fets_eval = fets_eval,
                       level = fe_level1 )
#        fe_grid1.deactivate( (1,0) )
#        fe_grid1.deactivate( (1,1) )
        
    fe_xdomain = XFESubDomain( domain = fe_domain,
                               fets_eval = xfets_eval,
                               #fe_grid_idx_slice = fe_grid1[1,0],
                                fe_grid_slice = fe_grid['X  -  0.5 '] )

    fe_xdomain.deactivate_sliced_elems()

    tstepper = TS( sdomain = fe_domain,
                   bcond_list =  [ BCDofGroup( var='u', value = 0., dims = [0],
                                               get_dof_method = fe_grid.get_left_dofs ),
#                                   BCDofGroup( var='u', value = 0., dims = [1],
#                                  get_dof_method = fe_grid.get_bottom_dofs ),                                  
                         BCDofGroup( var='u', value = .005, dims = [0],
                                  get_dof_method = fe_grid.get_top_right_dofs ) ],

         rtrace_list =  [ 
#                     RTraceDomainListField(name = 'Damage' ,
#                                    var = 'omega', idx = 0,
#                                    record_on = 'update',
#                                    warp = True),
#                     RTraceDomainListField(name = 'Displacement' ,
#                                    var = 'u', idx = 0,
#                                    record_on = 'update',
#                                    warp = True),
                                    #                    RTraceDomainListField(name = 'N0' ,
#                                      var = 'N_mtx', idx = 0,
#                                      record_on = 'update')
                ]             
            )
    
    # Add the time-loop control
    #global tloop
    tloop = TLoop( tstepper = tstepper, debug = True,
                   tline = TLine( min = 0.0,  step = 1.0, max = 1.0 ) )

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


    #import cProfile
    #cProfile.run('tloop.eval()', 'tloop_prof' )
    print tloop.eval()
    #import pstats
    #p = pstats.Stats('tloop_prof')
    #p.strip_dirs()
    #print 'cumulative'
    #p.sort_stats('cumulative').print_stats(20)
    #print 'time'
    #p.sort_stats('time').print_stats(20)

    # Put the whole thing into the simulation-framework to map the
    # individual pieces of definition into the user interface.
    #
    from ibvpy.plugins.ibvpy_app import IBVPyApp
    app = IBVPyApp( ibv_resource = tloop )
    app.main()    
    
if __name__ == '__main__':
    example_with_new_domain()