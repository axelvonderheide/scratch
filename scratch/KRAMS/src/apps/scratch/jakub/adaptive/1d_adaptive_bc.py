'''
Created on Mar 21, 2011

@author: jakub
'''
from ibvpy.core.astrategy import AStrategyBase
from ibvpy.api import BCDof

class ChangeBC( AStrategyBase ):

    bc_list = [BCDof( var = 'u', dof = 0, value = 0. ),
               BCDof( var = 'u', dof = 2, value = 1. ),
               BCDof( var = 'u', dof = 1, value = 1. ),
               BCDof( var = 'u', dof = 1, value = 1. )# redundant will be excluded by the procedure
               ]

    def begin_time_step( self, t ):
        '''Prepare a new load step
        '''
        self.tloop.bcond_list = self.bc_list
        print 'bc list '
        print self.bc_list
        K = self.tstepper.K
        K.reset()
        for bc in self.tloop.bcond_list:
            print 'constraint'
            bc.apply_essential( K )

    def end_time_step( self, t ):
        '''Prepare a new load step
        '''
        self.tloop.bcond_list = []
        self.bc_list.pop()


def example_with_new_domain():
    from ibvpy.api import \
        TStepper as TS, RTraceGraph, RTraceDomainListField, TLoop, \
        TLine, BCDof, IBVPSolve as IS, DOTSEval
    from ibvpy.mats.mats1D.mats1D_elastic.mats1D_elastic import MATS1DElastic
    from ibvpy.fets.fets1D.fets1D2l import FETS1D2L

    fets_eval = FETS1D2L( mats_eval = MATS1DElastic( E = 10. ) )
    from ibvpy.mesh.fe_grid import FEGrid

    # Discretization
    domain = FEGrid( coord_max = ( 2., 0., 0. ),
                           shape = ( 2, ),
                           fets_eval = fets_eval )

    ts = TS( dof_resultants = True,
             sdomain = domain,
         # conversion to list (square brackets) is only necessary for slicing of 
         # single dofs, e.g "get_left_dofs()[0,1]"

         bcond_list = [
#                       BCDof( var = 'u', dof = 0, value = 0. ),
#                       BCDof( var = 'u', dof = 2, value = 1. )
                        ],
         rtrace_list = [RTraceDomainListField( name = 'Stress' ,
                         var = 'sig_app', idx = 0 ),
                         RTraceDomainListField( name = 'Displacement' ,
                                    var = 'u', idx = 0,
                                    warp = True ),
                ]
            )

    # Add the time-loop control
    tloop = TLoop( tstepper = ts, debug = True,
                   adap = ChangeBC(),
                   tline = TLine( min = 0.0, step = 1., max = 2.0 ) )

    print '---- result ----'
    print tloop.eval()

    # Put the whole stuff into the simulation-framework to map the
    # individual pieces of definition into the user interface.
    #
#    from ibvpy.plugins.ibvpy_app import IBVPyApp
#    app = IBVPyApp( ibv_resource = tloop )
#    app.main()


if __name__ == '__main__':
    example_with_new_domain()
