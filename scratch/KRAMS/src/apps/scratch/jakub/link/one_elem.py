'''
Created on Feb 12, 2010

@author: jakub
'''

def example_with_new_domain():
    from ibvpy.api import \
        TStepper as TS, RTraceGraph, RTraceDomainListField, \
        RTraceDomainListInteg, TLoop, \
        TLine, BCDof, IBVPSolve as IS, DOTSEval, BCSlice
    from ibvpy.fets.fets2D.fets2D4q import FETS2D4Q
    from ibvpy.mats.mats2D.mats2D_elastic.mats2D_elastic import MATS2DElastic
    from ibvpy.mats.mats2D.mats2D_sdamage.mats2D_sdamage import MATS2DScalarDamage

    from ibvpy.api import BCDofGroup
    fets_eval = FETS2D4Q(mats_eval = MATS2DElastic(E= 1. ,nu= 0.)) 
    #fets_eval = FETS2D4Q(mats_eval = MATS2DScalarDamage()) 

    print fets_eval.vtk_node_cell_data
    
    from ibvpy.mesh.fe_grid import FEGrid
    from ibvpy.mesh.fe_refinement_grid import FERefinementGrid
    from ibvpy.mesh.fe_domain import FEDomain
    from mathkit.mfn.mfn_line.mfn_line import MFnLineArray

    # Discretization
    fe_grid = FEGrid( coord_max = (1.,1.,0.), 
                      shape   = (1, 1),
                      fets_eval = fets_eval )

    tstepper = TS( sdomain = fe_grid,
                   bcond_list =  [ 
                                  BCSlice(var='u', value = 0., dims = [0,1],
                                        slice =  fe_grid[0,0,0,:]),
#                                   BCSlice(var='u', value = 0., dims = [0],
#                                        link_slice = fe_grid[0,0,-1,-1],
#                                        link_coeffs = [1.],
#                                        link_dims = [0],
#                                        slice =  fe_grid[0,0,-1,0]),
                                    BCDofGroup( var='u', value = 0., dims = [0],
                                                get_dof_method = fe_grid.get_top_right_dofs,
                                                link_coeffs = [1.],
                                                link_dims = [0],
                                                get_link_dof_method = fe_grid.get_bottom_right_dofs ),
                                    BCSlice(var='f', value = 1., dims = [0],
                                        slice =  fe_grid[0,0,-1,-1]),
                                   ],
         rtrace_list =  [ 
                     RTraceDomainListField(name = 'Displacement' ,
                                    var = 'u', idx = 0,
                                    record_on = 'update',
                                    warp = True),
                ]             
            )
    
    # Add the time-loop control
    #global tloop
    tloop = TLoop( tstepper = tstepper, KMAX = 300, tolerance = 1e-4, debug = True,
                   tline = TLine( min = 0.0,  step = 1.0, max = 1.0 ) )

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