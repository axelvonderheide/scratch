from ibvpy.core.api import \
    TStepper as TS, RTraceGraph,\
    TLoop, TLine, BCDof,\
    IBVPSolve as IS
    
from ibvpy.api import MGridDomain, RTraceDomainField
    
from ibvpy.mats.mats1D.mats1D_damage.mats_damage1d import MATS1DDamage

from ibvpy.fets.fets1D.fets1D2l import FETS1D2L
from averaging import UniformDomainAveraging, LinearAF, QuarticAF
    
# Tseval for a discretized line domain
#
tseval  = UniformDomainAveraging( fets_eval = FETS1D2L(mats_eval = MATS1DDamage(E = 1.,
                                                                            epsilon_0 = 1.,
                                                                            epsilon_f = 1.,
                                                                            #stiffness = "algorithmic")) )
                                                                            stiffness = "secant")) ,
                                    correction = True,
                                    avg_function = LinearAF(Radius = 0.44))
# Discretization
#

    
from ibvpy.mesh.mgrid_domain import MeshGridAdaptor

mgrid_adaptor = MeshGridAdaptor( n_nodal_dofs = 1,
                                 n_e_nodes_geo = (1,0,0),
                                 n_e_nodes_dof = (1,0,0),
                                 node_map_geo = [0,1],
                                 node_map_dof = [0,1]
                                 )

line_domain = MGridDomain( lengths = (1.,0.,0.),
                              shape = (5,0,0),
                              adaptor = mgrid_adaptor )

   
    # Put the tseval (time-stepper) into the spatial context of the
# discretization and specify the response tracers to evaluate there.
#
right_dof = 5
ts = TS( tse = tseval,
     sdomain = line_domain,
     bcond_list = [ BCDof(var='u', dof = 0, value = 0.),
                BCDof(var='u', dof = 5, value = 1. ),
                BCDof(var='f', dof = 1, value = 1.e-9 ),
                ],
     rtrace_list = [ 
                RTraceGraph(name = 'Fi,right over u_right (iteration)' ,
                  var_y = 'F_int', idx_y = right_dof,
                  var_x = 'U_k', idx_x = right_dof,
                  record_on = 'update'),
#                  RTraceDomainField(name = 'Strain Field' ,
#                  var = 'eps', idx = 0,
#                  record_on = 'update'),
#                  RTraceDomainField(name = 'Strain Field' ,
#                  var = 'omega', idx = 0,
#                  record_on = 'update'),
#                  RTraceElemField(name = 'Deformation' ,
#                               var = 'u', idx = 0,
#                               record_on = 'update'),                
#              RTraceDomainField(name = 'Displacement Field' ,
#              var = 'u', idx = 0,
#              record_on = 'update'),
              RTraceDomainField(name = 'Damage' ,
              var = 'omega', idx = 0,
              record_on = 'update'),
         ]             
     )

# Add the time-loop control
#
tl = TLoop( tstepper = ts,
     DT = 1.,
     tline  = TLine( min = 0.0,  max = 4. ),
     tolerance = 1e-8)

if __name__ == '__main__':
    tl.eval()
    # Put the whole stuff into the simulation-framework to map the
    # individual pieces of definition into the user interface.
    #
    from ibvpy.plugins.ibvpy_app import IBVPyApp
    app = IBVPyApp( tloop = tl )
    app.main()