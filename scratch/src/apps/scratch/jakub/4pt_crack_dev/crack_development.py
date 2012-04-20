'''
Created on Mar 28, 2011

@author: jakub
'''
if __name__ == '__main__':
    def example_2d():
        from ibvpy.core.astrategy import AStrategyBase
        from ibvpy.api import BCDof, BCDofGroup, FEGrid, BCSlice, FEDomain, \
            FERefinementGrid, TStepper as TS, RTraceDomainListField
        from ibvpy.fets.fets2D.fets2Dtf import FETS2DTF
        from ibvpy.fets.fets2D.fets2D4q import FETS2D4Q
        from ibvpy.fets.fets2D.fets2D4q8u import FETS2D4Q8U
        from ibvpy.fets.fets2D.fets2D4q9u import FETS2D4Q9U
        from ibvpy.fets.fets2D.fets2D4q12u import FETS2D4Q12U
        from ibvpy.fets.fets2D.fets2D4q16u import FETS2D4Q16U
        from apps.scratch.jakub.global_enrichement.fets_crack_tf import FETSCrackTF
        from ibvpy.mesh.xfe_subdomain import XFESubDomain
        from enthought.traits.api import HasTraits, Array
        from ibvpy.mesh.fe_ls_domain import FELSDomain

        from ibvpy.mats.mats2D5.mats2D5_bond.mats2D5_plastic_bond import MATS2D5PlasticBond

        from ibvpy.core.tloop import TLoop, TLine
        import numpy as np

        #########################################
        #Material Parameters
        #########################################
        #E_m = 34 000 MPa
        #nu_m = 0.2
        #E_f = 54 000 MPa
        #Tau_fr = 3 MPA
        #G = 30 000 MPa
        #s_crit  assumed 1.e-5 yielding 

        m_eval = MATS2D5PlasticBond( E_m = 340. , #MPa including thickness 0.01M
                                     nu_m = 0.2, #
                                     E_f = 6., #MPa EA
                                     nu_f = 0.,
                                     G = 1.25e5 , #MPa 
                                     #sigma_y = 1.25 , #MPa \tau_fr
                                     sigma_y = 1.25 * 1000000. , #temporary set higher to avoid nonlinear comp.
                                     stress_state = 'plane_stress' )

        fets_eval = FETS2DTF( parent_fets = FETS2D4Q16U(),
                              mats_eval = m_eval )
        xfets_eval = FETSCrackTF( parent_fets = fets_eval, int_order = 5 ,
                                  tri_subdivision = 1 )

        H = lambda x: np.array( 0.5 * np.sign( x ) + 1, dtype = int )
        # Discretization
        fineness_x = 19
        fineness_y = 5

        fe_domain = FEDomain()
        fe_level1 = FERefinementGrid( domain = fe_domain, fets_eval = fets_eval )

        fe_grid1 = FEGrid( 
                           coord_max = ( 0.285, .076, 0. ),
                           shape = ( fineness_x, fineness_y ),
                           fets_eval = fets_eval,
                           level = fe_level1 )

#        #
#        #first crack from right: Y < 0.051 in 1.LS; 1 in 2. LS
#        fe_xdomain = XFESubDomain( domain = fe_domain,
#                                    fets_eval = xfets_eval,
#                                    boundary = lambda X, Y: H( -Y + 0.051 ) * H( X - 0.225 ), #
#                                    slice = fe_grid1['19047.61905005*X**4 -20215.92603891*X**3 +\
#                                                                8022.94611889*X**2 -1409.39497482*X+\
#                                                                        92.40622701-Y'] )#polyfit p=4
#        fe_xdomain.deactivate_sliced_elems()


#        fe_xdomain = XFESubDomain( domain = fe_domain,
#                                    fets_eval = xfets_eval,
#                                    boundary = lambda X, Y: H( -Y + 0.051 ) * H( X - 0.225 ), #
#                                    slice = fe_grid1['19047.61905005*X**4 -20215.92603891*X**3 +\
#                                                                8022.94611889*X**2 -1409.39497482*X+\
#                                                                        92.40622701-Y'] )#polyfit p=4

        bls_function = lambda X, Y: ( np.sign( X - 0.22 ) + np.sign( -Y + 0.051 ) ) - 1
        ls_function = lambda X, Y:19047.61905005 * X ** 4 - 20215.92603891 * X ** 3 + \
                        8022.94611889 * X ** 2 - 1409.39497482 * X + \
                                92.40622701 - Y

#        ls_function = lambda X, Y: X - 0.151 - Y
#        bls_function = lambda X, Y:-( ( X - 0.151 ) ** 2 + ( Y - 0.0 ) ** 2 - 0.05 ** 2 )

        fe_xdomain = FELSDomain( domain = fe_domain,
                                 fets_eval = xfets_eval,
                                 fe_grid = fe_grid1,
                                 ls_function = ls_function,
                                 bls_function = bls_function, #
                                 )

        fe_tip_xdomain = FELSDomain( domain = fe_domain,
                                     fets_eval = xfets_eval,
                                     fe_grid = fe_xdomain,
                                     ls_function = bls_function,
                                     )

        # deactivation must be done only after the dof enumeration has been completed
        fe_xdomain.deactivate_intg_elems_in_parent()
        fe_tip_xdomain.deactivate_intg_elems_in_parent()

        print 'dof offset arr ', fe_domain.dof_offset_arr

        #third crack from right: Y < 0.046 in 2.LS; 1 in 3.LS; X > 0.15
#        fe_xdomain2 = XFESubDomain( domain = fe_domain, #third crack
#                                   fets_eval = xfets_eval,
#                                   fe_grid_slice = fe_grid1['-88.81571047 *X**3+  56.5274648 *X**2 -\
#                                                            10.97186991*X +0.66858495-Y'] )#p=3
#        fe_xdomain2.deactivate_sliced_elems()
        #
        #fifth crack: 1 in 3. LS; X > 0.0875
        #fe_xdomain3 = XFESubDomain( domain = fe_domain, #fifth crack
        #                           fets_eval = xfets_eval,
        #                           fe_grid_slice = fe_grid1['9.50347099e+02*X**4 -5.47856002e+02*X**3+\
        #                                                     1.13942734e+02*X**2  -9.37037987*X+\
        #                                                    2.56078567e-01-Y'] )#p=4
        #fe_xdomain3.deactivate_sliced_elems()
        ##
        ##second crack: 0.035 in 2. and further LS
        #fe_xdomain4 = XFESubDomain( domain = fe_domain, #second crack
        #                           fets_eval = xfets_eval,
        #                           fe_grid_slice = fe_grid1['2.30769231*X**2+  0.09807692*X -\
        #                                     0.12504808 - Y  '] ) #p=2
        #fe_xdomain4.deactivate_sliced_elems()
        ##
        ##fourth crack: 0.02 in 3. and further LS
        #fe_xdomain5 = XFESubDomain( domain = fe_domain, #fourth crack
        #                       fe_grid_slice = fe_grid1['X - 0.615*0.25  - Y  '] ) #p=1
        #fe_xdomain5.deactivate_sliced_elems()
        #
        ##
        ##sixth crack: 0.051 in 3.LS; 0.066 in 4.LS; X > 0.0425
        #fe_xdomain6 = XFESubDomain( domain = fe_domain, #sixth crack
        #                       fets_eval = xfets_eval,
        #                       fe_grid_slice = fe_grid1['7.36512162e+02*X**4 -3.12383247e+02*X**3 +\
        #                                     4.72881197e+01*X**2  -2.41064925*X+\
        #                                    3.73082528e-02 - Y'] ) #p=4 
        #fe_xdomain6.deactivate_sliced_elems()
        #
        ##
        ##seventh crack: 0.025 in 4.LS; X > 0.0125
        #fe_xdomain7 = XFESubDomain( domain = fe_domain, #seventh crack
        #                       fets_eval = xfets_eval,
        #                       fe_grid_slice = fe_grid1[' 1.89166054e+02*X**3 +  6.13523745*X**2+\
        #                                   0.107873027*X  -7.68782666e-03 -Y'] ) #p=3


#        class FakeSlice( HasTraits ):
#            dofs = Array( int )
#            elems = Array( int )
#            dof_X = Array( float )
#
#        fs1 = FakeSlice()
#        c1 = fe_xdomain.elem_dof_map[-7:, -32:].reshape( -1, 16, 2 )
#        #c1 = fe_xdomain.elem_dof_map[:, -32:].reshape( -1, 16, 2 )
#        #c2a = fe_xdomain2.elem_dof_map[:-11, -32:].reshape( -1, 16, 2 )
#
#        #fs1.dofs = np.vstack( ( c1, c2a ) )
#        fs1.dofs = c1
#        fs1.elems = np.zeros( fs1.dofs.shape[0] )
#        fs1.dof_X = np.zeros( ( fs1.dofs.shape[0], 3 ) )
##        #
##        #
##        fs2 = FakeSlice()
##        c2b = fe_xdomain2.elem_dof_map[-11:, -32:].reshape( -1, 16, 2 )
##        fs2.dofs = c2b
##        fs2.elems = np.zeros( fs2.dofs.shape[0] )
##        fs2.dof_X = np.zeros( ( fs2.dofs.shape[0], 3 ) )


        bc_list = [
                   BCSlice( var = 'u', value = 0.0004 , dims = [0, 2], #value = 0.002
                             slice = fe_grid1[ -1, :-1, : ] ),
                    BCSlice( var = 'u', value = 0., dims = [0, 1],
                             slice = fe_grid1[ 0, 0, 0, 0 ] ),
                    BCSlice( var = 'u', value = -.001 , dims = [1],
                            slice = fe_grid1[-1, :, -1, :] ),

                    #redundant bc removed automatically during initiation
                    #BCDof( var = 'u', value = 0., dof = 1 )
                    ]

        boundary_list = [bls_function, bls_function,
                         #redundant bc removed automatically during initiation
                         bls_function
                    ]

        class ChangeBC( AStrategyBase ):


            def begin_time_step( self, t ):
                '''Prepare a new load step
                '''
                print 'begin_time_step'
                ls = boundary_list.pop()
                fe_xdomain.bls_function = ls

                fe_tip_xdomain.ls_function = ls

                fe_xdomain.deactivate_intg_elems_in_parent()
                fe_tip_xdomain.deactivate_intg_elems_in_parent()

                cdofs = fe_tip_xdomain.elem_xdof_map.flatten()

                const_list = [ BCDof( var = 'u', dof = dof, value = 0.0 )
                              for                   dof in cdofs ]

                print 'dof offset arr ', fe_domain.dof_offset_arr

                self.tloop.bcond_list = bc_list + const_list
                K = self.tstepper.K
                K.reset()

                for bc in self.tloop.bcond_list:
                    bc.setup( 'sctx' )  #BCDofGroup and BCSlice are realized during their setup
                                        #hack: sctx not used, therefore dummy string parsed 

                for bc in self.tloop.bcond_list:
                    bc.apply_essential( K )

                ##treatment of ls
                #fe_xdomain.boundary = boundary_list[0]

            def end_time_step( self, t ):
                '''Prepare a new load step
                '''
                print 'end_time_step'
                self.tloop.bcond_list = []
                #bc_list.pop()

        cdofs = fe_tip_xdomain.elem_xdof_map.flatten()
        bc_tip_list = [ BCDof( var = 'u', dof = dof, value = 0.0 )
                       for                   dof in cdofs ]

        ts = TS( dof_resultants = True,
                 sdomain = fe_domain,
                 bcond_list = bc_list + bc_tip_list,
                 rtrace_list = [
#                                 RTraceGraph(name = 'Fi,right over u_right (iteration)' ,
#                                       var_y = 'F_int', idx_y = 0,
#                                       var_x = 'U_k', idx_x = 1),
#                            RTraceDomainListField(name = 'Stress' ,
#                                 var = 'sig_app', idx = 0, warp = True ),
############################################
                             RTraceDomainListField( name = 'Displ matrix' ,
                                            var = 'u_m' ,
                                            warp = True,
                                            warp_var = 'u_m'
                                            ),
#                            RTraceDomainListField( name = 'Displ reinf' ,
#                                                   warp = True,
#                                            warp_var = 'u_f',
#                                            var = 'u_f' ),
#                            RTraceDomainListField( name = 'Strain matrix' ,
#                                            var = 'eps_m', warp_var = 'u_m' ),
#                            RTraceDomainListField( name = 'Strain reinf' ,
#                                            var = 'eps_f', warp_var = 'u_f' ), #,
##                            RTraceDomainListField( name = 'Stress matrix' ,
##                                            var = 'sig_app_m' ),
##                            RTraceDomainListField( name = 'Stress reinf' ,
##                                            var = 'sig_app_f' ),
#                            RTraceDomainListField( name = 'slip' ,
#                                            var = 'eps_b', warp_var = 'u_f' ),
#                            RTraceDomainListField( name = 'bond stress' ,
#                                            var = 'sig_b', warp_var = 'u_f' ),
                        ]
                    )
#        
#        # Add the time-loop control
        tloop = TLoop( tstepper = ts,
                       #adap = ChangeBC(),
                       tline = TLine( min = 0.0, step = .25, max = .5 ) )

        tloop.eval()

        from ibvpy.plugins.ibvpy_app import IBVPyApp
        ibvpy_app = IBVPyApp( ibv_resource = ts )
        ibvpy_app.main()


    example_2d()
