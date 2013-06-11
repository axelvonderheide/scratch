#-------------------------------------------------------------------------------
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
        TLine, IBVPSolve as IS, DOTSEval, BCSlice, FEGrid, BCDofGroup
    from ibvpy.fets.fets1D5.fets1D52l4uLRH import \
        FETS1D52L4ULRH, MATS1DElastic, MATS1DPlastic, MATS1D5Bond
        
    from ibvpy.fets.fets2D.fets2Drotsym import FETS2Drotsym
    from ibvpy.fets.fets2D.fets2D4q import FETS2D4Q
    from ibvpy.mats.mats2D.mats2D_elastic.mats2D_elastic import MATS2DElastic
    from mathkit.mfn.mfn_line.mfn_line import \
        MFnLineArray
    from math import sqrt, pi as Pi

    # Concrete parameters
    E_concrete = 32000.0               # [MPa]  E-Modulus of concrete
    nu_concrete = 0.2                  # [-]    Poisson ratio

    # Yarn parameters
    A_epoxy  = 1.760                   # [mm^2] cross-sectional area of epoxy
    A_glass  = 0.896                   # [mm^2] cross-sectional area of glass
    A_yarn   = A_epoxy + A_glass       # [mm^2] cross-sectional area of yarn
    P_yarn   = sqrt( 4 * A_yarn * Pi ) # [mm] perimeter of the yarn
    E_yarn   = 17000.0                 # [MPa]  effective E-Modulus of the impregnated yarn
    stiffness_fiber   = E_yarn * A_yarn
    
    penalty_stiffness = 1.e6
    
    # Bond parameters
    tau_max  = 13.0             # [N/mm^2] - frictional shear stress
    T_max    = tau_max * P_yarn # [N/mm] - frictional shear flow
    s_crit   = 0.028            # [mm] - onset of inelastic slip
    G        = T_max / s_crit   # [N/mm] - shear flow stiffness
    
    # Geometry
    L_e      = 5.                            # [mm] - embedded length
    
    # Loading conditions
    u_max    = 0.3
    # f_max   = 1200
    
    # Discretization
    fineness_x = 4
    fineness_y = 4
    
    # Material model construction
    mats_eval_bond = MATS1D5Bond( mats_phase1 = MATS1DElastic( E = stiffness_fiber ),
                                  mats_phase2 = MATS1DElastic( E = 0 ),
                                  mats_ifslip = MATS1DElastic( E = G ),
                                  mats_ifopen = MATS1DElastic( E = penalty_stiffness))
    
    mats_eval_matrix = MATS2DElastic( E = E_concrete,
                                      nu = nu_concrete,
                                      stress_state="plane_strain")

    # Finite element construction
    fets_eval_bond = FETS1D52L4ULRH( mats_eval = mats_eval_bond )
    fets_eval_matrix = FETS2Drotsym(parent_fets = FETS2D4Q(),
                                    mats_eval = MATS2DElastic(E=2,nu= .2,
                                    stress_state= 'rotational_symetry')) 

    # Discretization
    domain_bond = FEGrid( coord_min = (0., -L_e/5.),
                          coord_max = (L_e,  0.),
                          shape   = (fineness_x,1),
                          fets_eval = fets_eval_bond )
    
    domain_matrix1 = FEGrid( coord_max = (L_e, L_e),
                             shape   = (fineness_x,fineness_y),
                             fets_eval = fets_eval_matrix )
    
#    domain_matrix2 = FEGrid( coord_min = (L_e, 0.),
#                             coord_max = (L_e*2, L_e),
#                             shape   = (fineness_x,fineness_y),
#                             fets_eval = fets_eval_matrix )

    end_dof = domain_bond[-1,0,-1,0 ].dofs[0,0,0]
    ts = TS( dof_resultants = True,
             sdomain = [ domain_bond, domain_matrix1],#, domain_matrix2 ],
             bcond_list =  [                                
                            # Fixed fiber at the left end
                            BCSlice( var='u', value = 0., dims = [0],
                                     slice = domain_bond[0,0,0,0] ),

                            # Fixed y-displacement - no crack opening
                            BCSlice( var='u', value = 0., dims = [1],
                                     slice = domain_bond[:,0,:,0] ),
                                 
                            # Loading at the right end of the fiber
                            BCSlice( var='u', value = u_max, dims = [0],
                                     slice = domain_bond[-1,0,-1,0] ), 
                        
                            # Support the matrix in the horizontal direction
                            BCSlice( var='u', value = 0., dims = [0],
                                     slice = domain_matrix1[0,:,0,:]),
                                 
#                            # Support the matrix in the vertical direction
#                            BCSlice( var='u', value = 0., dims = [1],
#                                     slice = domain_matrix1[:,0,:,0]),
                        
#                            # Connect bond and matrix domains
#                            BCSlice( var='u', value = 0., dims = [0,1],
#                                     link_slice = domain_bond[:,-1,:,-1],
#                                     slice = domain_matrix1[:,0,:,0],
#                                     link_coeffs = [1.]),

                            BCDofGroup( var='u', value = 0., dims = [0,1],
                                        get_link_dof_method = domain_bond.get_top_dofs,
                                        get_dof_method = domain_matrix1.get_bottom_dofs,
                                        link_coeffs = [1.] ),
                                               
#                            # Connect two matrix domains
#                            BCSlice( var='u', value = 0., dims = [0,1],
#                                     slice = domain_matrix1[-1,:,-1,:],
#                                     link_slice = domain_matrix2[0,:,0,:],
#                                     link_coeffs = [1.]),
                                    
                        ],
         rtrace_list =  [ RTraceGraph(name = 'Fi,right over u_right (iteration)' ,
                               var_y = 'F_int', idx_y = end_dof,
                               var_x = 'U_k', idx_x = end_dof),
                          RTraceDomainListField(name = 'slip' ,
                                var = 'slip', idx = 0 ),
                          RTraceDomainListField(name = 'eps1' ,
                                var = 'eps1', idx = 0 ),
                          RTraceDomainListField(name = 'eps2' ,
                                var = 'eps2', idx = 0 ),
                          RTraceDomainListField(name = 'shear_flow' ,
                                var = 'shear_flow', idx = 0 ),
                          RTraceDomainListField(name = 'sig1' ,
                                var = 'sig1', idx = 0 ),
                          RTraceDomainListField(name = 'sig2' ,
                                var = 'sig2', idx = 0 ),
                          RTraceDomainListField(name = 'Displacement' ,
                                var = 'u', idx = 0), 
                          RTraceDomainListField(name = 'Stress' ,
                                var = 'sig_app', idx = 0)                   
                ] )
    
    # Add the time-loop control
    tloop = TLoop( tstepper = ts, KMAX = 30, debug = False,
                   tline  = TLine( min = 0.0,  step = 1.0, max = 1.0 ))
    
    print tloop.eval()    
    # Put the whole stuff into the simulation-framework to map the
    # individual pieces of definition into the user interface.
    from ibvpy.plugins.ibvpy_app import IBVPyApp
    app = IBVPyApp( ibv_resource = tloop )
    app.main()
if __name__ == '__main__':
    run()