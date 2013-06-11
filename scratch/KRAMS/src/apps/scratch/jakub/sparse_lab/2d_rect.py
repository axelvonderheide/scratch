

#----------------------- example --------------------
if __name__ == '__main__':
    from core.api import \
        TStepper as TS, MGridDomain, RTraceGraph, RTraceDomainField,  \
        TLine, BCDof, IBVPSolve as IS
        
    from sparse_dots_eval import SpDOTSEval
    from sparse_tloop import STLoop
    #from core.api import TLoop, DOTSEval
    #from lib.mats.mats2D.mats2D_cmdm.mats_mdm2d import MACMDM
    from lib.mats.mats2D.mats2D_elastic.mats2D_elastic import MATS2DElastic
    from lib.fets.fets2D.fets2D4q9u import FETS2D4Q9U

    fets_eval = FETS2D4Q9U(mats_eval = MATS2DElastic())            

    # Tseval for a discretized line domain
    #
    tseval  = SpDOTSEval( fets_eval = fets_eval )
        
    # Discretization
    #
    domain = MGridDomain( lengths = (1.,1.,0.), 
                             shape = (1,1,0), 
                             n_nodal_dofs = 2, 
                             # NOTE: the following properties must be defined and 
                             # must correspond to the used element formulation
                             n_e_nodes_geo = (1,1,0), 
                             n_e_nodes_dof = (2,2,0), 
                             node_map_geo = [0,1,3,2], 
                             node_map_dof = [0,2,8,6,1,5,7,3,4] )
                                 
        

#        for i, elem in enumerate( domain.elements ):
#            print 'elem.id_number', elem.id_number
#            print 'elem.nodes_geo', elem.nodes_geo
#            print 'elem.get_X_mtx()', elem.get_X_mtx()        
        
        # Put the tseval (time-stepper) into the spatial context of the
        # discretization and specify the response tracers to evaluate there.
        #
        
    right_dof = 2
    ts = TS( tse = tseval,
         sdomain = domain,
         # conversion to list (square brackets) is only necessary for slicing of 
         # single dofs, e.g "get_left_dofs()[0,1]"
         bcond_list =  [ BCDof(var='u', dof = i, value = 0.) for i in  domain.get_left_dofs()[:,0]  ] +
                    [ BCDof(var='u', dof = i, value = 0.) for i in [domain.get_left_dofs()[0,1]] ] +    
                    [ BCDof(var='u', dof = i, value = 0.002 ) for i in domain.get_right_dofs()[:,0] ],
         rtrace_list =  [ 
#                        RTraceGraph(name = 'Fi,right over u_right (iteration)' ,
#                               var_y = 'F_int', idx_y = right_dof,
#                               var_x = 'U_k', idx_x = right_dof,
#                               record_on = 'update'),
#                         RTraceDomainField(name = 'Stress' ,
#                         var = 'sig_app', idx = 0,
#                         record_on = 'update'),
                     RTraceDomainField(name = 'Displacement' ,
                                    var = 'u', idx = 0,
                                    record_on = 'update',
                                    warp = True),
#                             RTraceDomainField(name = 'N0' ,
#                                          var = 'N_mtx', idx = 0,
#                                          record_on = 'update')
                      
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
    #sim.configure_traits()
