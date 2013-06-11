'''
1D Bar
'''
if __name__ == '__main__':
    from core.api import \
        TStepper as TS, MGridDomain, RTraceGraph, RTraceDomainField,\
        TLine, BCDof,\
        IBVPSolve as IS, DOTSEval
        
    from sparse_tloop import STLoop
        
    #from lib.mats.mats1D.mats1D_damage.mats_damage1d import MATS1DDamage
    from lib.mats.mats1D.mats1D_elastic.mats1D_elastic import MATS1DElastic
    from lib.fets.fets1D.fets1D2l import FETS1D2L
        
    # Tseval for a discretized line domain
    #
    
    fets = FETS1D2L(mats_eval = MATS1DElastic())
    tseval  = DOTSEval( fets_eval = fets)
    # Discretization
    #

    line_domain = MGridDomain( lengths = (1., 0., 0.),
                                  shape = (1, 0, 0),
                                  n_e_nodes_geo = (1, 0, 0), 
                                  n_e_nodes_dof = (1, 0, 0), 
                                  node_map_geo = [0, 1], 
                                  node_map_dof = [0, 1],
                                   )

   
    # Put the tseval (time-stepper) into the spatial context of the
    # discretization and specify the response tracers to evaluate there.
    #
    right_dof = 1
    ts = TS( tse = tseval,
         sdomain = line_domain,
         bcond_list = [ BCDof(var='u', dof = 0, value = 0.),
                    BCDof(var='u', dof = 1, value = 1. ),
                    #BCDof(var='f', dof = 1, value = 1.e-6 ),
                    ],
         rtrace_list = [ 
                    RTraceGraph(name = 'Fi,right over u_right (iteration)' ,
                  var_y = 'F_int', idx_y = right_dof,
                  var_x = 'U_k', idx_x = right_dof,
                  record_on = 'update'),
#                  RTraceDomainField(name = 'Strain Field' ,
#                  var = 'eps', idx = 0,
#                  record_on = 'update'),
                  RTraceDomainField(name = 'Displacement Field' ,
                  var = 'u', idx = 0,
                  record_on = 'update'),
#                  RTraceDomainField(name = 'Damage Field' ,
#                  var = 'omega', idx = 0,
#                  record_on = 'update'),           
             ]             
         )

    # Add the time-loop control
    #
    tl = TLoop( tstepper = ts,
         DT = .001,
         #T  = TLine( min = 0.0,  max = 1. ),
         tolerance = 1e-8)
    
    tl.eval()
    # Put the whole stuff into the simulation-framework to map the
    # individual pieces of definition into the user interface.
    #
    sim = IS( tloop = tl )
    sim.configure_traits()
    
