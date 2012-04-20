
from enthought.traits.api import \
     Array, Bool, Callable, Enum, Float, HasTraits, Interface, implements, \
     Instance, Int, Trait, Str, Enum, Callable, List, TraitDict, Any, \
     on_trait_change, Tuple, WeakRef, Delegate, Property, cached_property, Dict

from enthought.traits.ui.api import \
     Item, View

from enthought.traits.ui.menu import \
     OKButton, CancelButton
     

from numpy import \
     zeros, float_, ix_, meshgrid

from ibvpy.api import \
     ITStepperEval, TStepperEval, DOTSEval, DOTSListEval

from ibvpy.core.rtrace_eval import RTraceEval
from ibvpy.fets.fets_eval import IFETSEval
from mathkit.matrix_la.sys_mtx_assembly import SysMtxAssembly


from time import time

if __name__ == '__main__':

    from ibvpy.api import \
        TStepper as TS, RTraceGraph, RTraceDomainField, TLoop, \
        TLine, BCDof, IBVPSolve as IS, DOTSEval, BCDofGroup
    from ibvpy.mats.mats1D.mats1D_elastic.mats1D_elastic import MATS1DElastic
    from ibvpy.mesh.fe_grid import FEGrid
    from ibvpy.mesh.fe_domain_list import FEDomainList
    from ibvpy.fets.fets1D.fets1D2l import FETS1D2L
    
    fets_eval = FETS1D2L(mats_eval = MATS1DElastic(E=10., A=1.))        

    # Discretization
    fe_domain1 = FEGrid( coord_max = (3.,0.,0.), 
                               shape   = (3,),
                               n_nodal_dofs = 1,
                               dof_r = fets_eval.dof_r,
                               geo_r = fets_eval.geo_r )

    fe_domain2 = FEGrid( coord_min = (3.,0.,0.),  
                               coord_max = (6.,0.,0.), 
                               shape   = (3,),
                               n_nodal_dofs = 1,
                               dof_r = fets_eval.dof_r,
                               geo_r = fets_eval.geo_r )

    fe_domain  = FEDomainList( subdomains = [ fe_domain1, fe_domain2 ] )
    #print "n_dofs " ,fe_domain.n_dofs
    # Tseval for a discretized line domain
    ts_eval1  = DOTSEval( fets_eval = fets_eval )
    ts_eval2  = DOTSEval( fets_eval = fets_eval )

    ts_eval = DOTSListEval( dots_list = [ts_eval1, ts_eval2] )
    
#    print "one ",fe_domain1.get_left_dofs()
#    print "two ",fe_domain2.get_left_dofs()
#    print "three ",fe_domain2.get_right_dofs()
#    
    ts = TS( tse = ts_eval,
             dof_resultants = True,
             sdomain = fe_domain,
             bcond_list = [
                            BCDofGroup( var='u', value = 0., dims = [0],
                                  get_dof_method = fe_domain1.get_left_dofs ),
                            BCDofGroup( var='u', value = 0., dims = [0],
                                  get_dof_method = fe_domain2.get_left_dofs,
                                    #link_dofs = [3],
                                    get_link_dof_method = fe_domain1.get_right_dofs,
                                    link_coeffs = [1.]),                                  
                            BCDofGroup( var='f', value = 1., dims = [0],
                                  get_dof_method = fe_domain2.get_right_dofs ) ],
#                            [BCDof(var='u', dof = 0, value = 0.),
#                            BCDof(var='u', dof = 4, link_dofs = [3], link_coeffs = [1.],
#                                  value = 0. ),
#                            BCDof(var='f', dof = 7, value = 1) ],
             rtrace_list =  [ RTraceGraph(name = 'Fi,right over u_right (iteration)' ,
                                   var_y = 'F_int', idx_y = 0,
                                   var_x = 'U_k', idx_x = 1),
#                        RTraceDomainField(name = 'Stress' ,
#                             var = 'sig_app', idx = 0),
#                         RTraceDomainField(name = 'Displacement' ,
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
    print ts.F_int
    print ts.rtrace_list[0].trace.ydata
    