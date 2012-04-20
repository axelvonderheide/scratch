'''
Created on Mar 22, 2011

@author: jakub
'''
from ibvpy.core.astrategy import AStrategyBase
from ibvpy.api import BCDof, BCDofGroup, FEGrid, BCSlice

from ibvpy.mats.mats2D.mats2D_elastic.mats2D_elastic import MATS2DElastic
from ibvpy.fets.fets2D.fets2D4q import FETS2D4Q
from ibvpy.fets.fets2D.fets2D4q8u import FETS2D4Q8U
from ibvpy.fets.fets2D.fets2D4q9u import FETS2D4Q9U
from ibvpy.fets.fets2D.fets2D9q import FETS2D9Q


fets_eval = FETS2D4Q( mats_eval = MATS2DElastic( E = 1., nu = 0. ) )


# Discretization
fe_grid1 = FEGrid( coord_max = ( 1., 1., 0. ),
                    shape = ( 1, 1 ),
                    fets_eval = fets_eval )

#print '------BC FORMAT---------'
#print fe_grid1[-1, 0, -1, 0].dofs


bc_list = [
#               BCDof( var = 'u', dof = 0, value = 0. ),
#               BCDof( var = 'u', dof = 1, value = 0. ),
#               BCDof( var = 'u', dof = 4, value = 1. ),
#               BCDof( var = 'u', dof = 5, value = 0. ),

#               BCDofGroup( var = 'u', dims = [0, 1], get_dof_method = fe_grid1.get_bottom_left_dofs , value = 0. ),
#               BCDofGroup( var = 'u', dims = [0], get_dof_method = fe_grid1.get_bottom_right_dofs, value = 1. ),
#               BCDofGroup( var = 'u', dims = [1], get_dof_method = fe_grid1.get_bottom_right_dofs, value = 0. ),

               BCSlice( var = 'u', dims = [0, 1], slice = fe_grid1[0, 0, 0, 0] , value = 0. ),
               BCSlice( var = 'u', dims = [1], slice = fe_grid1[-1, 0, -1, 0], value = 0. ),

               BCSlice( var = 'u', dims = [0], slice = fe_grid1[-1, 0, -1, 0], value = .1 ),
               BCSlice( var = 'u', dims = [0], slice = fe_grid1[-1, 0, -1, -1], value = .1 ),
               BCDof( var = 'u', dof = 1, value = 1. )# redundant will be excluded by the procedure
               ]


class ChangeBC( AStrategyBase ):

    def begin_time_step( self, t ):
        '''Prepare a new load step
        '''
        self.tloop.bcond_list = bc_list
        K = self.tstepper.K
        K.reset()
        for bc in self.tloop.bcond_list:
            bc.setup( 'sctx' )  #BCDofGroup and BCSlice are realized during their setup
                                #hack: sctx not used, therefore dummy string parsed 

        for bc in self.tloop.bcond_list:
            print 'constraint'
            bc.apply_essential( K )

    def end_time_step( self, t ):
        '''Prepare a new load step
        '''
        self.tloop.bcond_list = []
        bc_list.pop()


def example_with_new_domain():
    from ibvpy.api import \
        TStepper as TS, RTraceGraph, RTraceDomainListField, TLoop, \
        TLine, BCDof, IBVPSolve as IS, DOTSEval

#    fets_eval = FETS2D4Q( mats_eval = MATS2DElastic( E = 1., nu = 0. ) )
#
#
#    # Discretization
#    fe_grid1 = FEGrid( coord_max = ( 1., 1., 0. ),
#                    shape = ( 1, 1 ),
#                    fets_eval = fets_eval )


    ts = TS( dof_resultants = True,
             sdomain = fe_grid1,
         # conversion to list (square brackets) is only necessary for slicing of 
         # single dofs, e.g "get_left_dofs()[0,1]"

         bcond_list = [
                        ],
         rtrace_list = [
#                        RTraceDomainListField( name = 'Stress' ,
#                         var = 'sig_app', idx = 0 ),
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
    from ibvpy.plugins.ibvpy_app import IBVPyApp
    app = IBVPyApp( ibv_resource = tloop )
    app.main()


if __name__ == '__main__':
    example_with_new_domain()
