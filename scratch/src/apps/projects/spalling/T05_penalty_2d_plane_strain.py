#--------------------------------------------
# Copyright (c) 2009, IMB, RWTH Aachen.
# All rights reserved.
#
# This software is provided without warranty under the terms of the BSD
# license included in simvisage/LICENSE.txt and may be redistributed only
# under the conditions described in the aforementioned license.  The license
# is also available online at http://www.simvisage.com/licenses/BSD.txt
#
# Thanks for using Simvisage open source!
#
# Created on Sep 10, 2009 by: jakub
def run():
    from ibvpy.api import \
        TStepper as TS, RTraceGraph, RTraceDomainListField, TLoop, \
        TLine, IBVPSolve as IS, DOTSEval, BCSlice, FEGrid, BCDofGroup, \
        RTraceDomainListInteg
    from ibvpy.fets.fets1D5.fets1D52l4uLRH import \
        FETS1D52L4ULRH, MATS1DElastic, MATS1DPlastic, MATS1D5Bond
#    from ibvpy.fets.fets1D5.fets1D52l6uLRH import FETS1D52L6ULRH
#    from ibvpy.fets.fets1D5.fets1D52l8uLRH import FETS1D52L8ULRH    
    #from ibvpy.fets.fets1D.fets1D2l import FETS1D2L

    from ibvpy.fets.fets2D.fets2Drotsym import FETS2Drotsym
    from ibvpy.fets.fets2D.fets2D4q import FETS2D4Q
#    from ibvpy.fets.fets2D.fets2D4q8u import FETS2D4Q8U
#    from ibvpy.fets.fets2D.fets2D4q9u import FETS2D4Q9U
#    from ibvpy.fets.fets2D.fets2D4q12u import FETS2D4Q12U
    from ibvpy.mats.mats2D.mats2D_elastic.mats2D_elastic import MATS2DElastic
#    from ibvpy.mats.mats3D.mats3D_cmdm.mats3D_cmdm import \
#    MATS3DMicroplaneDamage, PhiFnStrainSoftening
    from ibvpy.mats.mats3D.mats3D_elastic.mats3D_elastic import MATS3DElastic
    from mathkit.mfn.mfn_line.mfn_line import \
        MFnLineArray
    from math import sqrt, pi as Pi

    # Concrete parameters
    E_concrete = 32000.0               # [MPa]  E-Modulus of concrete
    nu_concrete = 0.2                  # [-]    Poisson ratio

    # Yarn parameters
    A_epoxy = 1.760                   # [mm^2] cross-sectional area of epoxy
    A_glass = 0.896                   # [mm^2] cross-sectional area of glass
    A_yarn = A_epoxy + A_glass       # [mm^2] cross-sectional area of yarn
    P_yarn = sqrt( 4 * A_yarn * Pi ) # [mm] perimeter of the yarn
    E_yarn = 17000.0                 # [MPa]  effective E-Modulus of the impregnated yarn
    stiffness_fiber = E_yarn * A_yarn

    penalty_stiffness = 1.e6

    # Bond parameters
    #tau_max  = 13.0             # [N/mm^2] - frictional shear stress
    tau_max = 5.0
    T_max = tau_max * P_yarn # [N/mm] - frictional shear flow
    s_crit = 0.009                # [mm] - onset of inelastic slip
    G = T_max / s_crit    # [N/mm] - shear flow stiffness

    # Geometrynu_concrete
    L_e = 5.                       # [mm] embedded length

    # Loading conditions
    u_max = 0.005         # [mm]       
    # f_max   = 1200 

    # Discretization
    fineness_x = 4
    fineness_y = 4

    # Material model construction
    mats_eval_bond = MATS1D5Bond( mats_phase1 = MATS1DElastic( E = stiffness_fiber ),
                                  mats_phase2 = MATS1DElastic( E = 0. ),
                                  #mats_ifslip = MATS1DPlastic( E = G, sigma_y = T_max, K_bar = 0.,H_bar = 0.), # plastic function of slip
                                  mats_ifslip = MATS1DElastic( E = G ), # elastic function of slip
                                  #mats_ifopen = MATS1DElastic( E =penalty_stiffness)  # elastic function of open
                                  mats_ifopen = MATS1DElastic( stress_strain_curve = MFnLineArray( ydata = [ penalty_stiffness, 0., 0. ],
                                                                                                                                                      xdata = [ -1., 0., 1.0] ) )
                                   )
                                                                                                                                    # piecewise elastic function of open
#    mats_eval_matrix = MATS3DMicroplaneDamage( model_version = 'stiffness',
#                                   E = 32e3,
#                                   nu = 0.2,
#                                   phi_fn = PhiFnStrainSoftening( G_f = 42.e-3,  #N/mm,
#                                                                  f_t   = 4.0 )#N/mm^2
#                                   ) # the values are from Dissertation of Stephan Voss <Ingenieurmodelle zum Tragverhalten textilbewehrtem Beton> 

#    mats_eval_matrix = MATS3DElastic(E=1.e10, 
#                                                                nu=0.1)
    mats_eval_matrix = MATS2DElastic( E = E_concrete,
                                      nu = nu_concrete,
                                      stress_state = "plane_strain" )

    # Finite element construction
    fets_eval_bond = FETS1D52L4ULRH( mats_eval = mats_eval_bond )

#    fets_eval_matrix = FETS2Drotsym(
#                                    parent_fets = FETS2D4Q(),
#                                    mats_eval = mats_eval_matrix) 
    fets_eval_matrix = FETS2D4Q( mats_eval = mats_eval_matrix )

    # Discretization
    domain_bond = FEGrid( coord_min = ( 0., -L_e / 5. ),
                          coord_max = ( L_e, 0. ),
                          shape = ( fineness_x * 1, 1 ), # Caution! the bond and the matrix should match! Test each time!!!
                          fets_eval = fets_eval_bond )

    domain_matrix1 = FEGrid( 
                             coord_min = ( 0., 0. ),
                             coord_max = ( L_e, L_e ),
                             shape = ( fineness_x, fineness_y ),
                             fets_eval = fets_eval_matrix )

    end_dof = domain_bond[-1, 0, -1, 0 ].dofs[0, 0, 0]
    ts = TS( dof_resultants = True,
             sdomain = [ domain_bond, domain_matrix1],
             bcond_list = [
                            # Fixed fiber at the left end
                            BCSlice( var = 'u', value = 0., dims = [0],
                                     slice = domain_bond[0, 0, 0, 0] ),

                            # Fixed y-displacement - no crack opening
                            BCSlice( var = 'u', value = 0., dims = [1],
                                     slice = domain_bond[:, 0, :, 0] ),

                            # Support the matrix in the horizontal direction
                            BCSlice( var = 'u', value = 0., dims = [0],
                                     slice = domain_matrix1[0, :, 0, :] ),

                            BCSlice( var = 'u', value = 0., dims = [1],
                                     slice = domain_matrix1[0, 0, 0, 0] ),

                            # Loading at the right  top of the matrix
                            BCSlice( var = 'u', value = -u_max, dims = [1],
                                     slice = domain_matrix1[-1, -1, -1, -1] ),

#                            # Support the matrix in the vertical direction
#                            BCSlice( var='u', value = 0., dims =[1],
#                                     slice = domain_matrix1[:,0,:,0]),
#                        
#                            # Connect bond and matrix domains
                             BCDofGroup( var = 'u', value = 0., dims = [0, 1],
                                        get_link_dof_method = domain_bond.get_top_dofs,
                                        get_dof_method = domain_matrix1.get_bottom_dofs,
                                        link_coeffs = [1.] )

                        ],
         rtrace_list = [
                          RTraceGraph( name = 'Fi,right over u_right (iteration)' ,
                               var_y = 'F_int', idx_y = end_dof,
                               var_x = 'U_k', idx_x = end_dof ),
                          RTraceDomainListField( name = 'slip' ,
                                var = 'slip', idx = 0 ),
#                          RTraceDomainListField(name = 'eps1' ,
#                                var = 'eps1', idx = 0 ),
#                          RTraceDomainListField(name = 'eps2' ,
#                                var = 'eps2', idx = 0 ),
                          RTraceDomainListField( name = 'shear_flow' ,
                                var = 'shear_flow', idx = 0 ),
#                          RTraceDomainListField(name = 'sig1' ,
#                                var = 'sig1', idx = 0 ),
#                          RTraceDomainListField(name = 'sig2' ,
#                                var = 'sig2', idx = 0 ),
                          RTraceDomainListField( name = 'Displacement' ,
                                var = 'u', idx = 0, warp = True ),
#                          RTraceDomainListField(name = 'Stress' ,
#                                var = 'sig_app', idx = 0),   
                          RTraceDomainListField( name = 'fracture_energy' ,
                                        var = 'fracture_energy', idx = 0, warp = True,
                                        record_on = 'update' ),
#                          RTraceDomainListField(name = 'cohesive stress' ,
#                                var = 'cohesive_stress', idx = 0 ),
                          RTraceDomainListInteg( name = 'Total fracture energy' ,
                                        var = 'fracture_energy', idx = 0, warp = False,
                                        record_on = 'update' )
                ] )

    # Add the time-loop control
    tloop = TLoop( tstepper = ts, KMAX = 40, debug = False,
                   tline = TLine( min = 0.0, step = 1.0, max = 1.0 ) )

    print tloop.eval()
    # Put the whole stuff into the simulation-framework to map the
    # individual pieces of definition into the user interface.
    from ibvpy.plugins.ibvpy_app import IBVPyApp
    app = IBVPyApp( ibv_resource = tloop )
    app.main()
if __name__ == '__main__':
    run()
