from ibvpy.api import \
    BCDof, RTraceDomainField, TLoop, TLine

from ibvpy.core.ibv_model import IBVModel
from mgrid_domain import MeshGridAdaptor, MGridDomain, MeshManager, MGridDomain_enr
from ibvpy.fets.fets2D.fets2D4q import FETS2D4Q
from ibvpy.mats.mats2D.mats2D_elastic.mats2D_elastic import MATS2DElastic
from fets2D4q_disc import FETS2D4Q_disc
from dots_eval_enr import DOTSEval, DOTSManager
from tstepper_enr import TStepper_enr

class Multiscale(IBVModel):

    def _tloop_default(self):

        print 'tloop_default'

        # Define a mesh domain adaptor as a cached property to 
        # be constracted on demand
        mgrid_adaptor = MeshGridAdaptor( n_nodal_dofs = 2,
                                         # NOTE: the following properties must be defined and 
                                         # must correspond to the used element formulation
                                         n_e_nodes_geo = (1,1,0), 
                                         n_e_nodes_dof = (1,1,0), 
                                         node_map_geo = [0,1,3,2], 
                                         node_map_dof = [0,1,3,2] )

        # Discretization
        domain = MGridDomain( lengths = (2.,3.,0.), 
                              shape = (1,3,0), 
                              adaptor = mgrid_adaptor )
        
        fine_domain = MGridDomain_enr( lengths = (2.,1.,0.), 
                                   shape = (1,1,0), 
                                   adaptor = mgrid_adaptor,
                                   master = False)
        
        
        mmanager = MeshManager(mgrid_list =[domain,fine_domain])
                                         
        rtrace_list =  [ 
                         RTraceDomainField(name = 'Stress' ,
                                           var = 'sig_app', idx = 0,
                                           record_on = 'update'),
                         RTraceDomainField(name = 'Displacement' ,
                                           var = 'u', idx = 0,
                                           record_on = 'update',
                                           warp = True)
                         ]

        tseval  = DOTSEval( fets_eval = FETS2D4Q(mats_eval = MATS2DElastic(stress_state = "plane_stress",
                                                                                   E = 40.,
                                                                                   nu = 0.2)),
                                                                                   mesh = domain )
        fine_tseval  = DOTSEval( fets_eval = FETS2D4Q(mats_eval = MATS2DElastic(stress_state = "plane_stress",
                                                                                   E = 40.,
                                                                                   nu = 0.2)), mesh = fine_domain )
        dmanager = DOTSManager(dots_list = [tseval,fine_tseval])

        ts = TStepper_enr( tse = dmanager,
                       sdomain = mmanager,
                       bcond_list = [ BCDof(var='u', dof = i, value = 0.) for i in  domain.get_left_dofs()[:,0]  ] +
                    [ BCDof(var='u', dof = i, value = 0.) for i in domain.get_left_dofs()[:,1] ] +    
                    [ BCDof(var='u', dof = i, value = 0.002 ) for i in domain.get_top_right_dofs()[:,1] ]+
                    [ BCDof(var='u', dof = i, value = -0.002 ) for i in domain.get_bottom_right_dofs()[:,1] ],
                       rtrace_list = rtrace_list
                       )

 
            
        # Put the time-stepper into the time-loop
        #
        tmax = 1.
        # tmax = 0.0006
        n_steps = 1
        #n_steps = 3
    
        tloop = TLoop( tstepper = ts,
                    KMAX = 100, RESETMAX = 0,
                    tline = TLine( min = 0.0,  step=tmax/n_steps, max = tmax ) )

        return tloop

if __name__ == '__main__':
    from ibvpy.plugins.ibvpy_app import IBVPyApp
    mscale = Multiscale()
    
    app = IBVPyApp( ibv_resource = mscale )
    app.main()
