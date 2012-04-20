

from ibvpy.mats.mats2D.mats2D_elastic.mats2D_elastic import MATS2DElastic
from ibvpy.fets.fets2D.fets2D4q import FETS2D4Q
from ibvpy.fets.fets2D.fets2D4q8u import FETS2D4Q8U
from ibvpy.fets.fets_ls.fets_crack import FETSCrack
from ibvpy.api import\
     BCDofGroup, TStepper as TS, TLoop, TLine, RTraceGraph
from ibvpy.rtrace.rt_domain_list_field import RTraceDomainListField
from ibvpy.mesh.fe_grid_domain import FEGridDomain
from ibvpy.mesh.fe_subgrid_domain import FESubGridDomain
from ibvpy.mesh.fe_domain_list import FEDomainList
from ibvpy.mesh.fe_domain_tree import FEDomainTree

def combined_fe2D4q_with_fe2D4q8u():

    fets_eval_4u = FETS2D4Q(mats_eval = MATS2DElastic(E= 1.,nu = 0.))
    fets_eval_8u = FETS2D4Q8U(mats_eval = MATS2DElastic())
    xfets_eval = FETSCrack(parent_fets = fets_eval_4u) # should be set automatically

    # Discretization
    fe_domain1 = FEGridDomain( coord_max = (2.,6.,0.), 
                               shape   = (1,3),
                               fets_eval = fets_eval_4u )

    fe_subdomain = FESubGridDomain( parent_domain = fe_domain1,
                                    #fets_eval = fets_eval_8u,
                                    fets_eval = fets_eval_4u,
                                    #fets_eval = xfets_eval,
                                    fine_cell_shape = (1,1) )
    
    fe_subdomain.refine_elem( (0,1) )
    elem = fe_subdomain.elements
    m_elem = fe_domain1.elements
    print "nodes ",elem[0]

    fe_domain  = FEDomainList( subdomains = [ fe_domain1 ] )

    ts = TS( dof_resultants = True,
             sdomain = fe_domain,
             bcond_list =  [BCDofGroup(var='u', value = 1., dims = [1],
                                       get_dof_method = fe_domain1.get_top_dofs ),
                            BCDofGroup(var='u', value = 0., dims = [1],
                                       get_dof_method = fe_domain1.get_bottom_dofs ),
                            BCDofGroup(var='u', value = 0., dims = [0],
                                       get_dof_method = fe_domain1.get_bottom_left_dofs ),
                                       ],
             rtrace_list =  [ 
#                             RTraceGraph(name = 'Fi,right over u_right (iteration)' ,
#                                   var_y = 'F_int', idx_y = 0,
#                                   var_x = 'U_k', idx_x = 1),
                        RTraceDomainListField(name = 'Stress' ,
                             var = 'sig_app', idx = 1, warp = True ),
                        RTraceDomainListField(name = 'Displ' ,
                             var = 'u', idx = 1, warp = True ),                             
                        
#                             RTraceDomainField(name = 'Displacement' ,
#                                        var = 'u', idx = 0),
#                                 RTraceDomainField(name = 'N0' ,
#                                              var = 'N_mtx', idx = 0,
#                                              record_on = 'update')
                          
                    ]             
                )
    
    # Add the time-loop control
    tloop = TLoop( tstepper = ts,
                   tline  = TLine( min = 0.0,  step = 1, max = 1.0 ))
    
    print tloop.eval()
    print "nodes after",elem[0]
    from ibvpy.plugins.ibvpy_app import IBVPyApp
    ibvpy_app = IBVPyApp( ibv_resource = tloop )
    ibvpy_app.main()
        
#    print ts.F_int
#    print ts.rtrace_list[0].trace.ydata

if __name__ == '__main__':
    combined_fe2D4q_with_fe2D4q8u()