'''
Created on Sep 9, 2009

@author: jakub
'''

def cube():
    from ibvpy.api import \
        TStepper as TS, RTraceGraph, RTraceDomainListField, TLoop, \
        TLine, BCDofGroup, IBVPSolve as IS, DOTSEval, RTraceDomainListInteg
        
    from ibvpy.mats.mats3D.mats3D_elastic.mats3D_elastic import MATS3DElastic
    from ibvpy.fets.fets3D.fets3D8h import FETS3D8H
    from ibvpy.mats.mats3D.mats3D_cmdm.mats3D_cmdm import \
        MATS3DMicroplaneDamage, PhiFnStrainSoftening

    mats = MATS3DMicroplaneDamage( model_version = 'stiffness',
                                   E = 32.e3, #N/mm^2
                                   nu = 0.2,   
                                   phi_fn = PhiFnStrainSoftening( G_f = 42.e-3, #N/mm
                                                                                        f_t = 4.0))  #N/mm^2
    
    fets_eval = FETS3D8H(mats_eval = mats)            
                                            
    from ibvpy.mesh.fe_grid import FEGrid
    
    # Discretization
    domain = FEGrid( coord_max = (1.,1.,1.), 
                           shape   = (2, 2, 2),
                           fets_eval = fets_eval)
                                 
    ts = TS(
            sdomain = domain,
             bcond_list = [BCDofGroup( var='u', value = 0., dims = [0],
                                  get_dof_method = domain.get_left_dofs ),
                        BCDofGroup( var='u', value = 0., dims = [1,2],
                                  get_dof_method = domain.get_bottom_left_dofs ),                                  
                        BCDofGroup( var='u', value = 1.e-4, dims = [0],
                                  get_dof_method = domain.get_right_dofs ) ],
             rtrace_list = [ 
#                        RTraceGraph(name = 'Fi,right over u_right (iteration)' ,
#                                  var_y = 'F_int', idx_y = right_dof,
#                                  var_x = 'U_k', idx_x = right_dof,
#                                  record_on = 'update'),
                        RTraceDomainListField(name = 'Deformation' ,
                                       var = 'eps_app', idx = 0,
                                       record_on = 'update'),
                        RTraceDomainListField(name = 'fracture_energy' ,
                                        var = 'fracture_energy', idx = 0, warp = True,
                                        record_on = 'update'), 
                        RTraceDomainListInteg(name = 'Total fracture energy' ,
                                        var = 'fracture_energy', idx = 0, warp = False,
                                        record_on = 'update'),
                         RTraceDomainListField(name = 'Displacement' ,
                                        var = 'u', idx = 0),
                        RTraceDomainListField(name = 'Stress' ,
                         var = 'sig_app', idx = 0, warp = True,
                         record_on = 'update'),
#                         RTraceDomainListField(name = 'Stress' ,
#                                        var = 'sig', idx = 0,
#                                        record_on = 'update'),
#                        RTraceDomainListField(name = 'N0' ,
#                                       var = 'N_mtx', idx = 0,
#                                       record_on = 'update')
                        ]             
            )
    
    # Add the time-loop control
    tloop = TLoop( tstepper = ts,tolerance = 1e-5,
         tline  = TLine( min = 0.0,  step = 1., max = 1.0 ))
  
    tloop.eval()    
    
    # Put the whole stuff into the simulation-framework to map the
    # individual pieces of definition into the user interface.
    #
    from ibvpy.plugins.ibvpy_app import IBVPyApp
    app = IBVPyApp( ibv_resource = tloop )
    app.main()
    
def cylinder():
    from ibvpy.api import \
        TStepper as TS, RTraceGraph, RTraceDomainListField, \
        RTraceDomainListInteg, TLoop, \
        TLine, BCDof, IBVPSolve as IS, DOTSEval
    from ibvpy.mats.mats2D.mats2D_elastic.mats2D_elastic import MATS2DElastic
    from ibvpy.mats.mats2D.mats2D_sdamage.mats2D_sdamage import MATS2DScalarDamage
    from ibvpy.fets.fets2D.fets2Drotsym import FETS2Drotsym
    from ibvpy.fets.fets2D.fets2D4q import FETS2D4Q
    from ibvpy.fets.fets2D.fets2D4q8u import FETS2D4Q8U
    
    from ibvpy.api import BCDofGroup
    from ibvpy.mats.mats3D.mats3D_cmdm.mats3D_cmdm import \
    MATS3DMicroplaneDamage, PhiFnStrainSoftening
    from math import sqrt, pi
    
#    mats =  MATS2DElastic(E=2,nu= .2,
#                          stress_state= 'rotational_symetry')
    mats = MATS3DMicroplaneDamage( model_version = 'stiffness',
                                   E = 34e3,
                                   nu = 0.2,
                                   phi_fn = PhiFnStrainSoftening( G_f = 0.001117,
                                                                                        f_t = 2.8968 ))
                                   
    fets_eval = FETS2Drotsym(parent_fets = FETS2D4Q(),
                             mats_eval = mats)

    fets_eval.vtk_r *= 0.9
    from ibvpy.mesh.fe_grid import FEGrid
    from ibvpy.mesh.fe_refinement_grid import FERefinementGrid
    from ibvpy.mesh.fe_domain import FEDomain
    from mathkit.mfn.mfn_line.mfn_line import MFnLineArray

    radius = sqrt(1./pi)
#    f_i = (radius/2.)*2*pi
#    f_o = (radius)*2*pi
#    print 'f ',f_i,' ', f_o
    # Discretization
    fe_grid = FEGrid( #coord_min = (0.,radius/2.,0.), 
                     coord_max = (1.,radius,0.), 
                      shape   = (2, 1),
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
                         BCDofGroup( var='u', value = 1.e-4, dims = [0],
                                  get_dof_method = fe_grid.get_right_dofs ) ],
                                  
         rtrace_list =  [
                         RTraceDomainListField(name = 'Stress' ,
                         var = 'sig_app', idx = 0, warp = True,
                         record_on = 'update'),
                        RTraceDomainListField(name = 'fracture_energy' ,
                                        var = 'fracture_energy', idx = 0, warp = True,
                                        record_on = 'update'), 
                        RTraceDomainListInteg(name = 'Total fracture energy' ,
                                        var = 'fracture_energy', idx = 0, warp = False,
                                        record_on = 'update'),                                      
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
    tloop = TLoop( tstepper = tstepper,tolerance = 1e-5,
                   tline = TLine( min = 0.0,  step = 1.0, max = 1.0 ) )

    tloop.eval()

    from ibvpy.plugins.ibvpy_app import IBVPyApp
    app = IBVPyApp( ibv_resource = tloop )
    app.main()
    
if __name__ == '__main__':
    cylinder()