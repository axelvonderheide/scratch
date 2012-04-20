from ibvpy.api import \
    TStepper as TS, MGridDomain, RTraceGraph, TLoop, \
    TLine, BCDof, IBVPSolve as IS, DOTSEval

from ibvpy.api import MGridDomain, RTraceDomainField
    
from ibvpy.mats.mats2D.mats2D_sdamage.mats2D_sdamage import MATS2DScalarDamage
from ibvpy.mats.mats2D.mats2D_sdamage.strain_norm2d import Euclidean
from ibvpy.fets.fets2D.fets2D4q9u import  FETS2D4Q9U
from averaging import UniformDomainAveraging, LinearAF, QuarticAF

fets_eval = FETS2D4Q9U(mats_eval = MATS2DScalarDamage(E = 1.,nu = 0., 
                                               epsilon_0 = 1.e-3, 
                                               epsilon_f = 1.e-2,
                                               strain_norm_type = 'Euclidean'))            


from numpy import array, cos, sin, pi,sqrt
# Tseval for a discretized line domain
#
tseval  = UniformDomainAveraging( fets_eval = fets_eval,
                                 correction = True,
                                 avg_function = QuarticAF(Radius = 0.2)  )
    
# Discretization
#
from ibvpy.mesh.mgrid_domain import MeshGridAdaptor

mgrid_adaptor = MeshGridAdaptor( n_nodal_dofs = 2, 
                                 n_e_nodes_geo = (1,1,0), 
                                 n_e_nodes_dof = (2,2,0), 
                                 node_map_geo = [0,1,3,2], 
                                 node_map_dof = [0,2,8,6,1,5,7,3,4] )  

domain = MGridDomain( lengths = (2., 1., 0.), 
                         shape = (10, 5, 0), 
                         # NOTE: the following properties must be defined and 
                         # must correspond to the used element formulation
                         adaptor = mgrid_adaptor ) 

    

#        for i, elem in enumerate( domain.elements ):
#            print 'elem.id_number', elem.id_number
#            print 'elem.nodes_geo', elem.nodes_geo
#            print 'elem.get_X_mtx()', elem.get_X_mtx()        
    
# Put the tseval (time-stepper) into the spatial context of the
# discretization and specify the response tracers to evaluate there.
#
domain.n_dofs    

angle = 2.
angle_r = angle/180. * pi
s_angle = sin(angle_r)
c_angle = cos(angle_r)
diag = sqrt(2.)

ts = TS(
        tse = tseval,
        sdomain = domain,
         # conversion to list (square brackets) is only necessary for slicing of 
         # single dofs, e.g "get_left_dofs()[0,1]" which elsewise retuns an integer only
         bcond_list =  [ # constraint for all left dofs in x-direction: 
                     BCDof(var='u', dof = domain.get_bottom_left_dofs()[0,0], value = 0.)] +
                      # constraint for all left dofs in y-direction:
                    [BCDof(var='u', dof = domain.get_bottom_left_dofs()[0,1], value = 0.)] + 
                    [BCDof(var='u', dof = domain.get_bottom_right_dofs()[0,0], value = 2.*c_angle - 2.)]+
                    [BCDof(var='u', dof = domain.get_bottom_right_dofs()[0,1], value = 2.*s_angle)]+
                    [BCDof(var='u', dof = domain.get_top_left_dofs()[0,0], value = -1.*s_angle)]+ 
                    [BCDof(var='u', dof = domain.get_top_left_dofs()[0,1], value = c_angle - 1.)]+ 
                    [BCDof(var='u', dof = domain.get_top_right_dofs()[0,0], value = 2.*c_angle - s_angle - 2.)]+
                    [BCDof(var='u', dof = domain.get_top_right_dofs()[0,1], value = 2.*s_angle + c_angle - 1.)]
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
                                       var = 'u', idx = 1,
                                       record_on = 'update',
                                       warp = True),
                        RTraceDomainField(name = 'Damage' ,
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
         DT = .5,
         tline  = TLine( min = 0.0,  max = 1.0 ))
if __name__ == '__main__':  
    tl.eval()    
    # Put the whole stuff into the simulation-framework to map the
    # individual pieces of definition into the user interface.
    #
    from ibvpy.plugins.ibvpy_app import IBVPyApp
    app = IBVPyApp( tloop = tl )
    app.main()