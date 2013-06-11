'''
Created on Sep 8, 2009

@author: jakub
'''

def example_with_new_domain():
    from ibvpy.api import \
        TStepper as TS, RTraceGraph, RTraceDomainListField, \
        RTraceDomainListInteg, TLoop, \
        TLine, BCDof, IBVPSolve as IS, DOTSEval
    from ibvpy.mats.mats2D.mats2D_elastic.mats2D_elastic import MATS2DElastic
    from ibvpy.mats.mats2D.mats2D_sdamage.mats2D_sdamage import MATS2DScalarDamage
    from ibvpy.fets.fets2D.fets2D4q import FETS2D4Q
    from ibvpy.fets.fets2D.fets2D4q8u import FETS2D4Q8U
    from ibvpy.api import BCDofGroup
    from ibvpy.fets.fets2D.fets2Drotsym import FETS2Drotsym
    from ibvpy.mats.mats3D.mats3D_cmdm.mats3D_cmdm import \
    MATS3DMicroplaneDamage, PhiFnStrainSoftening
    me = MATS3DMicroplaneDamage( model_version = 'stiffness',
                                   E = 34e3,
                                   nu = 0.2,
                                   phi_fn = PhiFnStrainSoftening( G_f = 0.001117,
                                                                  f_t = 2.8968 ))
                                   
#    me = MATS2DElastic(E=2,nu= .2,
#                 stress_state= 'rotational_symetry')
    fets_eval = FETS2Drotsym(parent_fets = FETS2D4Q8U(),
                             mats_eval = me) 
    #fets_eval = FETS2D4Q(mats_eval = MATS2DScalarDamage()) 
    fets_eval.parent_fets.vtk_r *= 0.9
    
    from ibvpy.mesh.fe_grid import FEGrid
    from ibvpy.mesh.fe_refinement_grid import FERefinementGrid
    from ibvpy.mesh.fe_domain import FEDomain
    from mathkit.mfn.mfn_line.mfn_line import MFnLineArray

 
    # Discretization
    fe_grid = FEGrid( #coord_min = (radius/2.,0.,0.), 
                     coord_max = (0.01,0.01,0.), 
                      shape   = (2, 2),
                      fets_eval = fets_eval )

    

    tstepper = TS( sdomain = fe_grid,
                   bcond_list =  [ BCDofGroup( var='u', value = 0., dims = [0],
                                               get_dof_method = fe_grid.get_left_dofs ),
                                    BCDofGroup( var='u', value = 0., dims = [1],
                                               get_dof_method = fe_grid.get_bottom_left_dofs ),      
#                                   BCDofGroup( var='u', value = 0., dims = [1],
#                                  get_dof_method = fe_grid.get_bottom_dofs ), 
#                         BCDofGroup( var='u', value = 1., dims = [1],
#                                  get_dof_method = fe_grid.get_top_left_dofs ) ,                                 
                         BCDofGroup( var='u', value = .001, dims = [0],
                                  get_dof_method = fe_grid.get_top_right_dofs ) ],
                                  
         rtrace_list =  [
                         RTraceDomainListField(name = 'Stress' ,
                         var = 'sig_app', idx = 0, warp = True,
                         record_on = 'update'),
#                     RTraceDomainListField(name = 'Damage' ,
#                                    var = 'omega', idx = 0,
#                                    record_on = 'update',
#                                    warp = True),
                     RTraceDomainListField(name = 'Displacement' ,
                                    var = 'u', idx = 0,
                                    record_on = 'update',
                                    warp = True),
                                    #                    RTraceDomainListField(name = 'N0' ,
#                                      var = 'N_mtx', idx = 0,
#                                      record_on = 'update')
                ]             
            )
    
    # Add the time-loop control
    #global tloop
    tloop = TLoop( tstepper = tstepper, KMAX = 300, tolerance = 1e-4,
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
