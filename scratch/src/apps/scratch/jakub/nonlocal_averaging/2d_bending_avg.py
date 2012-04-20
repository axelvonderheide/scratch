'''
Created on Jan 27, 2010
        
@author: jakub
'''
from ibvpy.api import \
    TStepper as TS, RTraceGraph, TLoop, \
    TLine, BCDof, IBVPSolve as IS, DOTSEval, FEDomain, FERefinementGrid,\
    FEGrid, BCSlice


from apps.scratch.jakub.mlab.mlab_trace import RTraceDomainListField
from ibvpy.mats.mats2D.mats2D_sdamage.mats2D_sdamage import MATS2DScalarDamage
from ibvpy.mats.mats2D.mats2D_elastic.mats2D_elastic import MATS2DElastic
from ibvpy.mats.mats2D.mats2D_sdamage.strain_norm2d import Euclidean, Mazars, Rankine
from ibvpy.fets.fets2D.fets2D4q import FETS2D4Q
from ibvpy.fets.fets2D.fets2D4q9u import  FETS2D4Q9U
from ibvpy.fets.fets2D.fets2D4q8u import FETS2D4Q8U
from averaging import UniformDomainAveraging, LinearAF, QuarticAF
from numpy import array, cos, sin, pi,sqrt, deg2rad, arctan
from mathkit.mfn.mfn_line.mfn_line import MFnLineArray
from ibvpy.dots.avg_fn import AveragingFunction, LinearAF,QuarticAF

def app():
    mp = MATS2DScalarDamage(E = 34.e3,
                            nu = 0.2,
                            epsilon_0 = 1.e-2,
                            epsilon_f = 1.e-1,
                            #stiffness  = "algorithmic",
                            stiffness  = "secant",
                            strain_norm = Rankine())
    me = MATS2DElastic(E = 34e3,
                       nu = 0.2,
                       stress_state = "plane_strain")
        
    fets_eval = FETS2D4Q(mats_eval = mp)#, ngp_r = 3, ngp_s = 3)                                               
    # Discretization
        
    fe_domain = FEDomain()
    fe_level1 = FERefinementGrid( domain = fe_domain, 
                                  fets_eval = fets_eval,
                                  averaging = QuarticAF(radius = 0.3))
    fe_grid = FEGrid( coord_max = (2.,1.,0.), 
                      shape   = (20,10),
                      fets_eval = fets_eval,
                      level = fe_level1 )
            
    mf = MFnLineArray( #xdata = arange(10),
                       ydata = array([0,2,2.2]) )     
            

    
    ts = TS(sdomain = fe_domain,
             # conversion to list (square brackets) is only necessary for slicing of 
             # single dofs, e.g "get_left_dofs()[0,1]" which elsewise retuns an integer only
             bcond_list =  [
                        # constraint for all left dofs in y-direction:
                        BCSlice(var='u', slice = fe_grid[0,0,0,0],dims=[0,1], value = 0.), 
                        BCSlice(var='u', slice = fe_grid[-1,0,-1,0],dims=[0,1], 
                                time_function = mf.get_value, value = 0.),
                        BCSlice(var='u', slice = fe_grid[9,-1,-1,-1],dims=[1], 
                                time_function = mf.get_value, value = -0.02),
                        ],
             rtrace_list = [ 
    #                        RTraceGraph(name = 'Fi,right over u_right (iteration)' ,
    #                                  var_y = 'F_int', idx_y = right_dof,
    #                                  var_x = 'U_k', idx_x = right_dof,
    #                                  record_on = 'update'),
                            RTraceDomainListField(name = 'Deformation' ,
                                           var = 'eps', idx = 0,
                                           record_on = 'update'),
                            RTraceDomainListField(name = 'Displacement' ,
                                           var = 'u', idx = 1,
                                           record_on = 'update',
                                           warp = True),
                            RTraceDomainListField(name = 'Damage' ,
                                           var = 'omega', idx = 0,
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
                tolerance = 1.e-4,
                tline  = TLine( min = 0.0, step = 1., max = 7.0 ))
    tl.eval()    
    # Put the whole stuff into the simulation-framework to map the
    # individual pieces of definition into the user interface.
    #
    from ibvpy.plugins.ibvpy_app import IBVPyApp
    ibvpy_app = IBVPyApp( ibv_resource = ts )
    ibvpy_app.main()

if __name__ == '__main__':  
    app()