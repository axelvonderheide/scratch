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
        from mathkit.matrix_la.sys_mtx_assembly import SysMtxAssembly
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
                                     sigma_y = 1.0 * 50. , #temporary set higher to avoid nonlinear comp.
                                     stress_state = 'plane_stress' )

        fets_eval = FETS2DTF( parent_fets = FETS2D4Q(),
                              mats_eval = m_eval )
        xfets_eval = FETSCrackTF( parent_fets = fets_eval, int_order = 5 ,
                                  tri_subdivision = 0 )

        # Discretization
        fineness_x = 3
        fineness_y = 3

        fe_domain = FEDomain()
        fe_level1 = FERefinementGrid( domain = fe_domain, fets_eval = fets_eval )

        fe_grid1 = FEGrid( 
                           coord_max = ( 0.2, 0.2, 0. ),
                           shape = ( fineness_x, fineness_y ),
                           fets_eval = fets_eval,
                           level = fe_level1 )

        # inclined crack with moving circles as domain 
        ls_function = lambda X, Y: X - 0.07 - Y
        bls_function1 = lambda X, Y:-( ( X - 0.07 ) ** 2 + ( Y - 0.0 ) ** 2 - 0.08 ** 2 )
        bls_function2 = lambda X, Y:-( ( X - 0.07 ) ** 2 + ( Y - 0.0 ) ** 2 - 0.14 ** 2 )

        # vertical crack with lower halfspace as a domain
        ls_function = lambda X, Y: X - 0.0343 - Y
        ls_function2 = ls_function # lambda X, Y: X - 0.03 - Y
        bls_function1 = lambda X, Y:-( Y - 0.1734 )
        bls_function2 = lambda X, Y:-( Y - 0.09596 )
        bls_function3 = lambda X, Y:-( Y - 0.19 )

        fe_xdomain = FELSDomain( domain = fe_domain,
                                 fets_eval = xfets_eval,
                                 fe_grid = fe_grid1,
                                 ls_function = ls_function,
                                 bls_function = bls_function1, #
                                 )

        fe_tip_xdomain = FELSDomain( domain = fe_domain,
                                     fets_eval = xfets_eval,
                                     fe_grid = fe_xdomain,
                                     ls_function = bls_function1,
                                     )

        # deactivation must be done only after the dof enumeration has been completed
        fe_xdomain.deactivate_intg_elems_in_parent()
        fe_tip_xdomain.deactivate_intg_elems_in_parent()

        bc_list = [
#                    BCSlice( var = 'u', value = 0.0004 , dims = [0, 2], #value = 0.002
#                             slice = fe_grid1[ -1, :-1, : ] ),
                    BCSlice( var = 'u', value = 0., dims = [0, 1],
                             slice = fe_grid1[ 0, :, 0, : ] ),
                    BCSlice( var = 'u', value = .01 , dims = [0],
                            slice = fe_grid1[-1, :, -1, :] ),

                    #redundant bc removed automatically during initiation
                    #BCDof( var = 'u', value = 0., dof = 1 )
                    ]

        ls_list = [ls_function2, ls_function,
                         #redundant bc removed automatically during initiation
                         ls_function
                    ]

        boundary_list = [bls_function2, bls_function1,
                         #redundant bc removed automatically during initiation
                         bls_function1
                    ]

        class ChangeBC( AStrategyBase ):


            def begin_time_step( self, t ):
                '''Prepare a new load step
                '''

                # before changing anything - remember the old dof_offsets of the domains
                old_dof_offset_arr = fe_domain.dof_offset_arr.copy()

                print 'begin_time_step'
                bls = boundary_list.pop()
                ls = ls_list.pop()

                # reactivate everything
                fe_grid1.inactive_elems = []
                fe_xdomain.idx_mask[...] = False
                fe_tip_xdomain.idx_mask[...] = False

                fe_xdomain.ls_function = ls
                fe_xdomain.bls_function = bls
                fe_tip_xdomain.ls_function = bls

                fe_xdomain.deactivate_intg_elems_in_parent()
                fe_tip_xdomain.deactivate_intg_elems_in_parent()

                old_U_k = self.tloop.U_k.copy()

                self.tstepper.sdomain.changed_structure = True

                # Initialize the variables
                #
                self.tloop.R_k = self.tstepper.new_resp_var()
                self.tloop.U_k = self.tstepper.U_k
                self.tloop.d_U = self.tstepper.d_U
                self.tloop.U_n = self.tstepper.new_cntl_var()
                self.tstepper.K = SysMtxAssembly()

                new_U_k = self.tloop.U_k

                new_dof_offset_arr = fe_domain.dof_offset_arr
                #copy the segments of the old displacements into the new vector
                for os, oe, ns, ne in zip( old_dof_offset_arr[:-1], old_dof_offset_arr[1:],
                                           new_dof_offset_arr[:-1], new_dof_offset_arr[1:] ):
                    min_length = min( oe - os, ne - ns )
                    new_U_k[ ns: ns + min_length ] = old_U_k[ os: os + min_length ]

                self.cdofs = fe_tip_xdomain.elem_xdof_map.flatten()

                # set the crack tip (jump that was there so far to zero)
                #
                new_U_k[ self.cdofs ] = 0.0

                const_list = [ BCDof( var = 'u', dof = dof, value = 0.0 )
                              for                   dof in self.cdofs ]

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

        #cdofs = fe_tip_xdomain.elem_xdof_map.flatten()
        #bc_tip_list = [ BCDof( var = 'u', dof = dof, value = 0.0 )
        #               for                   dof in cdofs ]

        ts = TS( dof_resultants = True,
                 sdomain = fe_domain,
                 #bcond_list = bc_list + bc_tip_list,
                 rtrace_list = [
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
                            RTraceDomainListField( name = 'bond stress' ,
                                            var = 'sig_b', warp_var = 'u_f' ),
                        ]
                    )
#        
#        # Add the time-loop control
        adap = ChangeBC()
        tloop = TLoop( tstepper = ts,
                       adap = adap,
                       debug = False,
                       tline = TLine( min = 0.0, step = .25, max = .25 ) )

        u = tloop.eval()

        from ibvpy.plugins.ibvpy_app import IBVPyApp
        ibvpy_app = IBVPyApp( ibv_resource = ts )
        ibvpy_app.main()


    example_2d()
