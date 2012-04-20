from ibvpy.core.api import \
    TStepper as TS, RTraceGraph,  TLoop, \
    TLine, BCDof, IBVPSolve as IS
  
from ibvpy.api import MGridDomain, RTraceDomainField

from ibvpy.dots.dots_eval import DOTSEval
    
from ibvpy.mats.mats3D.mats3D_sdamage.mats3D_sdamage import MATS3DScalarDamage
from ibvpy.mats.mats3D.mats3D_sdamage.strain_norm3d import Euclidean
from ibvpy.fets.fets3D.fets3D8h import FETS3D8H

fets_eval = FETS3D8H(mats_eval = MATS3DScalarDamage(E = 1.,nu = 0., 
                                               epsilon_0 = 1.e-3, 
                                               epsilon_f = 1.e-2,
                                               strain_norm = Euclidean()))            
from averaging import UniformDomainAveraging, LinearAF

from numpy import array, cos, sin, pi,sqrt
# Tseval for a discretized line domain
#
tseval  = UniformDomainAveraging( fets_eval = fets_eval,
                                 correction = False,
                                 avg_function = LinearAF(Radius = 0.35)  )


from ibvpy.mesh.mgrid_domain import MeshGridAdaptor

# Define a mesh domain adaptor as a cached property to 
# be constracted on demand
mgrid_adaptor = MeshGridAdaptor( n_nodal_dofs = 3, 
                                 n_e_nodes_geo = (1, 1, 1), 
                                 n_e_nodes_dof = (1, 1, 1), 
                                 node_map_geo = [0, 1, 2, 3, 4, 5, 6, 7], 
                                 node_map_dof = [0, 1, 2, 3, 4, 5, 6, 7] )         
# Discretization
#
domain = MGridDomain( lengths = (1., 1., 1.), 
                         shape = (5, 5, 5), 
                         adaptor = mgrid_adaptor)
                             
    

#        for i, elem in enumerate( domain.elements ):
#            print 'elem.id_number', elem.id_number
#            print 'elem.nodes_geo', elem.nodes_geo
#            print 'elem.get_X_mtx()', elem.get_X_mtx()        
    
# Put the tseval (time-stepper) into the spatial context of the
# discretization and specify the response tracers to evaluate there.
#
domain.n_dofs    
cd_bot_fro = [domain.enum_nodes_dof[0, 0, 0], domain.enum_nodes_dof[-1, 0, 0]]
cd_top_fro = [domain.enum_nodes_dof[0, -1, 0], domain.enum_nodes_dof[-1, -1, 0]]
cd_top_bac = [domain.enum_nodes_dof[0, -1, -1], domain.enum_nodes_dof[-1, -1, -1]]
cd_bot_bac = [domain.enum_nodes_dof[0, 0, -1], domain.enum_nodes_dof[-1, 0, -1]]


dl_bot_fro = array([domain.nodes_dof[i].dofs  for i in cd_bot_fro])
dl_top_fro = array([domain.nodes_dof[i].dofs  for i in cd_top_fro])
dl_top_bac = array([domain.nodes_dof[i].dofs  for i in cd_top_bac])
dl_bot_bac = array([domain.nodes_dof[i].dofs  for i in cd_bot_bac])

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
                     BCDof(var='u', dof = i, value = 0.)\
                      for i in  dl_bot_fro[:,2]  ] +
                      # constraint for all left dofs in y-direction:
                    [ BCDof(var='u', dof = i, value = 0.)\
                      for i in  dl_bot_fro[:,1]] + 
                    [ BCDof(var='u', dof = dl_bot_fro[0,0], value = 0.)]+ 
                      # imposed displacement for all right dofs in x-direction:
                    [ BCDof(var='u', dof = i, value =  s_angle)\
                      for i in  dl_bot_bac[:,1]]+
                    [ BCDof(var='u', dof = i, value = (c_angle - 1.) )\
                      for i in  dl_bot_bac[:,2]]+
                    [ BCDof(var='u', dof = i, value = (s_angle+c_angle - 1))\
                      for i in  dl_top_bac[:,1]]+
                    [ BCDof(var='u', dof = i, value = (c_angle-s_angle - 1))\
                      for i in  dl_top_bac[:,2]]+
                    [ BCDof(var='u', dof = i, value =  (c_angle - 1))\
                      for i in  dl_top_fro[:,1]]+
                    [ BCDof(var='u', dof = i, value = -s_angle )\
                      for i in  dl_top_fro[:,2]] 
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