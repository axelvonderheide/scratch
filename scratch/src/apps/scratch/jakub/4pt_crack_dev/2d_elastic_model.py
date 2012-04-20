'''
Created on Mar 31, 2011

@author: jakub
'''
if __name__ == '__main__':
    def example_2d():
        from ibvpy.core.astrategy import AStrategyBase
        from ibvpy.api import BCDof, BCDofGroup, FEGrid, BCSlice, FEDomain, FERefinementGrid, \
            TStepper as TS, RTraceDomainListField
        from ibvpy.core.tloop import TLoop, TLine
        from ibvpy.fets.fets2D.fets2Dtf import FETS2DTF
        from ibvpy.fets.fets2D.fets2D4q import FETS2D4Q
        from ibvpy.fets.fets2D.fets2D4q8u import FETS2D4Q8U
        from ibvpy.fets.fets2D.fets2D4q9u import FETS2D4Q9U
        from ibvpy.fets.fets2D.fets2D4q12u import FETS2D4Q12U
        from ibvpy.fets.fets2D.fets2D4q16u import FETS2D4Q16U
        from ibvpy.rtrace.rt_dof import RTraceGraph

        from ibvpy.mats.mats2D5.mats2D5_bond.mats2D_bond import MATS2D5Bond
        from ibvpy.mats.mats2D5.mats2D5_bond.mats2D5_plastic_bond import MATS2D5PlasticBond
        from ibvpy.mats.mats2D.mats2D_elastic.mats2D_elastic import MATS2DElastic

        import numpy as np
        #########################################
        #Material Parameters
        #########################################
        #E_m = 34 000 MPa
        #nu_m = 0.2
        #E_f = 54 000 MPa
        #Tau_fr = 1.25 MPA
        #G = 12 500 MPa
        # s_crit  assumed 1.e-4 yielding 

        #########################################
        #Geometry
        #########################################
        #WEB
        ####
        #b = 18 mm
        #h = 84 mm 
        #
        ########
        #FLANGES
        ########
        #b = 110 mm
        #h = 18 mm (including the height of the transfer zone 6 mm)
        #

#        m_eval_tf = MATS2D5PlasticBond( E_m = 612. , #MPa including thickness 0.01M
#                                     nu_m = 0.2, #
#                                     E_f = 6., #MPa EA
#                                     nu_f = 0.,
#                                     G = 1.25e4 , #MPa 
#                                     #sigma_y = 1.25 , #MPa \tau_fr
#                                     sigma_y = 1.25 ,
#                                     stress_state = 'plane_stress' )

        m_eval_tf = MATS2D5Bond( E_m = 34000 * 0.018 , #MPa including thickness 0.01M
                                     nu_m = 0.2, #
                                     E_f = 54000. * 2.22e-4, #MPa EA
                                     nu_f = 0.,
                                     G = 1.25e4 , #MPa 
                                     stress_state = 'plane_stress' )


        m_eval_web = MATS2DElastic( E = 34000 * 0.018 + 54000. * 2.22e-4, #MPa including thickness 0.01M
                                    nu = 0.2, #
                                    stress_state = 'plane_stress' )

        m_eval_flange = MATS2DElastic( E = 34000 * 0.11 + 54000. * 1.854e-3, #MPa including thickness 0.01M
                                       nu = 0.2, #
                                       stress_state = 'plane_stress' )

        fets_eval_tf = FETS2DTF( parent_fets = FETS2D4Q8U(),
                                 mats_eval = m_eval_tf )

        fets_eval_web = FETS2D4Q8U( mats_eval = m_eval_web )

        fets_eval_flange = FETS2D4Q8U( mats_eval = m_eval_flange )


        # Discretization
        fineness_x_left = 4
        fineness_x_right = 2
        fineness_y_web = 1
        fineness_y_flange = 1

        fe_domain = FEDomain()
        fe_level1 = FERefinementGrid( domain = fe_domain, fets_eval = fets_eval_tf )
        ##
        #web
        ##

        fe_grid1 = FEGrid( coord_max = ( 0.285, .084 ),
                           shape = ( fineness_x_left, fineness_y_web ),
                           fets_eval = fets_eval_tf,
                           level = fe_level1
                           )
        fe_level2 = FERefinementGrid( domain = fe_domain, fets_eval = fets_eval_web )
        fe_grid2 = FEGrid( coord_min = ( 0.285, 0. ),
                           coord_max = ( 0.450, .084 ),
                           shape = ( fineness_x_right, fineness_y_web ),
                           fets_eval = fets_eval_web,
                           level = fe_level2
                           )

        ##
        #flanges
        ##
        ##lower
        fe_level3 = FERefinementGrid( domain = fe_domain, fets_eval = fets_eval_flange )
        fe_grid3 = FEGrid( coord_min = ( 0., -0.018 ),
                           coord_max = ( 0.285, 0. ),
                           shape = ( fineness_x_left, fineness_y_flange ),
                           fets_eval = fets_eval_flange,
                           level = fe_level3
                           )
        fe_level4 = FERefinementGrid( domain = fe_domain, fets_eval = fets_eval_flange )
        fe_grid4 = FEGrid( coord_min = ( 0.285, -0.018 ),
                           coord_max = ( 0.450, 0., ),
                           shape = ( fineness_x_right, fineness_y_flange ),
                           fets_eval = fets_eval_flange,
                           level = fe_level4
                           )
        ##upper
        fe_level5 = FERefinementGrid( domain = fe_domain, fets_eval = fets_eval_flange )
        fe_grid5 = FEGrid( coord_min = ( 0., 0.084 ),
                           coord_max = ( 0.285, 0.102 ),
                           shape = ( fineness_x_left, fineness_y_flange ),
                           fets_eval = fets_eval_flange,
                           level = fe_level5
                           )
        fe_level6 = FERefinementGrid( domain = fe_domain, fets_eval = fets_eval_flange )
        fe_grid6 = FEGrid( coord_min = ( 0.285, 0.084, ),
                           coord_max = ( 0.450, 0.102, ),
                           shape = ( fineness_x_right, fineness_y_flange ),
                           fets_eval = fets_eval_flange,
                           level = fe_level6
                           )

        def get_g2_bottom_dofs():
            a = np.sort( np.unique( np.hstack( ( fe_grid2[1:-1, 0, :, 0].dofs.flatten(),
                                               fe_grid2[0, 0, 1:-1, 0].dofs.flatten(),
                                               fe_grid2[-1, 0, :-1, 0].dofs.flatten() ) ) ) ).reshape( -1, 2 )
            return a, []

        def get_g2_top_dofs():
            a = np.sort( np.unique( np.hstack( ( fe_grid2[1:-1, -1, :, -1].dofs.flatten(),
                                               fe_grid2[0, -1, 1:, -1].dofs.flatten(),
                                               fe_grid2[-1, -1, :-1, -1].dofs.flatten() ) ) ) ).reshape( -1, 2 )
            return a, []

        def get_g4_top_dofs():
            a = np.sort( np.unique( np.hstack( ( fe_grid4[1:-1, -1, :, -1].dofs.flatten(),
                                               fe_grid4[0, -1, 1:, -1].dofs.flatten(),
                                               fe_grid4[-1, -1, :-1, -1].dofs.flatten() ) ) ) ).reshape( -1, 2 )
            return a, []

        def get_g6_bottom_dofs():
            a = np.sort( np.unique( np.hstack( ( fe_grid6[1:-1, 0, :, 0].dofs.flatten(),
                                               fe_grid6[0, 0, 1:-1, 0].dofs.flatten(),
                                               fe_grid6[-1, 0, :-1, 0].dofs.flatten() ) ) ) ).reshape( -1, 2 )
            return a, []


        bottom_right_dof = fe_grid4[-1, 0, -1, 0].dofs[0, 0, 1]

        ts = TS( #dof_resultants = True,
                 #sdomain = fe_domain,
                 sdomain = [fe_grid1, fe_grid2, fe_grid3, fe_grid4, fe_grid5, fe_grid6],
                 bcond_list = [
                               #RHS symmetry constraints 
                               BCDofGroup( var = 'u', value = 0.0 , dims = [0],
                                    get_dof_method = fe_grid2.get_right_dofs ),
                               BCDofGroup( var = 'u', value = 0.0 , dims = [0],
                                    get_dof_method = fe_grid4.get_right_dofs ),
                               BCDofGroup( var = 'u', value = 0.0 , dims = [0],
                                    get_dof_method = fe_grid6.get_right_dofs ),

                               # support 
                               BCDofGroup( var = 'u', value = 0., dims = [0, 1],
                                    get_dof_method = fe_grid3.get_bottom_left_dofs ),

                               #stitch together the layers of two field zone
                               BCDofGroup( var = 'u', value = 0., dims = [2, 3], link_dims = [0, 1],
                                    link_coeffs = [1.],
                                    get_link_dof_method = fe_grid1.get_right_dofs,
                                    get_dof_method = fe_grid1.get_right_dofs ),
                               BCDofGroup( var = 'u', value = 0., dims = [2, 3], link_dims = [0, 1],
                                    link_coeffs = [1.],
                                    get_link_dof_method = fe_grid1.get_top_dofs,
                                    get_dof_method = fe_grid1.get_top_dofs ),
                               BCDofGroup( var = 'u', value = 0., dims = [2, 3], link_dims = [0, 1],
                                    link_coeffs = [1.],
                                    get_link_dof_method = fe_grid1.get_bottom_dofs,
                                    get_dof_method = fe_grid1.get_bottom_dofs ),


                                # bottom flange -> web left
                                BCDofGroup( var = 'u', value = 0., dims = [0, 1], link_dims = [0, 1],
                                    link_coeffs = [1.],
                                    get_link_dof_method = fe_grid3.get_top_dofs,
                                    get_dof_method = fe_grid1.get_bottom_dofs ),

                                # web -> top flange left
                                BCDofGroup( var = 'u', value = 0., dims = [0, 1], link_dims = [0, 1],
                                    link_coeffs = [1.],
                                    get_dof_method = fe_grid1.get_top_dofs,
                                    get_link_dof_method = fe_grid5.get_bottom_dofs ),

                               #web cracked area-elastic area
                               BCDofGroup( var = 'u', value = 0., dims = [0, 1], link_dims = [2, 3],
                                    link_coeffs = [1.],
                                    get_link_dof_method = fe_grid1.get_right_dofs,
                                    get_dof_method = fe_grid2.get_left_dofs ),

                                #flange connection
                                BCDofGroup( var = 'u', value = 0., dims = [0, 1], link_dims = [0, 1],
                                    link_coeffs = [1.],
                                    get_dof_method = fe_grid3.get_right_dofs,
                                    get_link_dof_method = fe_grid4.get_left_dofs ),
                                BCDofGroup( var = 'u', value = 0., dims = [0, 1], link_dims = [0, 1],
                                    link_coeffs = [1.],
                                    get_dof_method = fe_grid5.get_right_dofs,
                                    get_link_dof_method = fe_grid6.get_left_dofs ),

                                #web -> bottom flange right
                                BCDofGroup( var = 'u', value = 0., dims = [0, 1], link_dims = [0, 1],
                                    link_coeffs = [1.],
                                    get_dof_method = get_g2_bottom_dofs,
                                    get_link_dof_method = get_g4_top_dofs ),

                               #web -> top flange right
                                BCDofGroup( var = 'u', value = 0., dims = [0, 1], link_dims = [0, 1],
                                    link_coeffs = [1.],
                                    get_dof_method = get_g2_top_dofs,
                                    get_link_dof_method = get_g6_bottom_dofs ),

                               #load
                               BCSlice( var = 'u', value = -.001 , dims = [1],
                                    slice = fe_grid6[-1, :, -1, :] ),
                               BCSlice( var = 'u', value = -.001 , dims = [1],
                                    slice = fe_grid2[-1, :, -1, :] ),
                               BCSlice( var = 'u', value = -.001 , dims = [1],
                                    slice = fe_grid4[-1, :, -1, :] ),
                               ],
                 rtrace_list = [#todo: create tracer integrating the forces along the right line
                                 RTraceGraph( name = 'Fi,right over u_right' ,
                                       var_y = 'F_int', idx_y = bottom_right_dof,
                                       var_x = 'U_k', idx_x = bottom_right_dof ),
                            RTraceDomainListField( name = 'Stress' ,
                                 var = 'sig_app', idx = 0, warp = True ),
                RTraceDomainListField( name = 'Displ matrix' ,
                                            var = 'u_m' ,
                                            warp = True,
                                            warp_var = 'u_m'
                                            ),
                RTraceDomainListField( name = 'Stress' ,
                                            var = 'sig_app' ),
                        ]
                    )
#        
#        # Add the time-loop control
        tloop = TLoop( tstepper = ts,
                       tline = TLine( min = 0.0, step = .25, max = .25 ) )


        tloop.eval()
        from ibvpy.plugins.ibvpy_app import IBVPyApp
        ibvpy_app = IBVPyApp( ibv_resource = ts )
        ibvpy_app.main()


    example_2d()
