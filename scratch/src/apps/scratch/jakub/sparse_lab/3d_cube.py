#----------------------- example --------------------

if __name__ == '__main__':
    from core.api import \
        TStepper as TS, MGridDomain, RTraceGraph, RTraceDomainField, \
        TLine, BCDof, IBVPSolve as IS, DOTSEval
        
    from sparse_tloop import STLoop

    from lib.mats.mats3D.mats3D_elastic.mats3D_elastic import MATS3DElastic
    from lib.fets.fets3D.fets3D8h import FETS3D8H
    from math import e

#    fets_eval = FETS2D9Q(mats_eval = MA2DSca    larDamage(strain_norm = Euclidean())) 
    fets_eval = FETS3D8H(mats_eval = MATS3DElastic())            

    # Tseval for a discretized line domain
    #
    tseval  = DOTSEval( fets_eval = fets_eval )
        
    # Discretization
    #
    domain = MGridDomain( lengths = (1., 1., 1.), 
                             shape = (10, 10, 10), 
                             # NOTE: the following properties must be defined and 
                             # must correspond to the used element formulation
                             n_nodal_dofs = 3, 
                             n_e_nodes_geo = (1, 1, 1), 
                             n_e_nodes_dof = (1, 1, 1), 
                             node_map_geo = [0, 1, 2, 3, 4, 5, 6, 7], 
                             node_map_dof = [0, 1, 2, 3, 4, 5, 6, 7] )
    
    domain.changed = True
                                 
        

#        for i, elem in enumerate( domain.elements ):
#            print 'elem.id_number', elem.id_number
#            print 'elem.nodes_geo', elem.nodes_geo
#            print 'elem.get_X_mtx()', elem.get_X_mtx()        
        
    # Put the tseval (time-stepper) into the spatial context of the
    # discretization and specify the response tracers to evaluate there.
    #
        
    right_dof = 2
    ts = TS(
            tse = tseval,
            sdomain = domain,
             # conversion to list (square brackets) is only necessary for slicing of 
             # single dofs, e.g "get_left_dofs()[0,1]" which elsewise retuns an integer only
             bcond_list =  [ # constraint for all left dofs in x-direction: 
                         BCDof(var='u', dof = i, value = 0.)\
                          for i in  domain.get_left_dofs()[:,0]  ] +
                          # constraint for all left dofs in y-direction:
                        [ BCDof(var='u', dof = i, value = 0.)\
                          for i in domain.get_left_dofs()[:,1]  ] + 
                          # imposed displacement for all right dofs in x-direction:
                        [ BCDof(var='u', dof = i, value = 2.*e-5 )\
                          for i in  domain.get_right_dofs()[:,0]]
                        ,
             rtrace_list = [ 
#                        RTraceGraph(name = 'Fi,right over u_right (iteration)' ,
#                                  var_y = 'F_int', idx_y = right_dof,
#                                  var_x = 'U_k', idx_x = right_dof,
#                                  record_on = 'update'),
#                        RTraceDomainField(name = 'Deformation' ,
#                                       var = 'eps', idx = 0,
#                                       record_on = 'update'),
                         RTraceDomainField(name = 'Displacement' ,
                                        var = 'u', idx = 0,
                                        record_on = 'update',
                                        warp = True),
#                         RTraceDomainField(name = 'Stress' ,
#                                        var = 'sig', idx = 0,
#                                        record_on = 'update'),
#                        RTraceDomainField(name = 'N0' ,
#                                       var = 'N_mtx', idx = 0,
#                                       record_on = 'update')
                        ]             
            )
    
    # Add the time-loop control
    #
    tl = TLoop( tstepper = ts,
         DT = 0.5)
         #T  = TLine( min = 0.0,  max = 1.0 ))
  
    tl.eval()    
    # Put the whole stuff into the simulation-framework to map the
    # individual pieces of definition into the user interface.
    #
    sim = IS( tloop = tl )
    sim.configure_traits()


