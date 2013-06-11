#-------------------------------------------------------------------------------
#
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
# Created on Sep 10, 2009 by: rch

def run():
    '''
    Pull-Out test with epoxy-impregnated yarn. Publisched in 
    
    Konrad, M., Chudoba, R., Tensile Behavior of Cementitous Composite Reinforced
    with Epoxy Impregnated Multifilament Yarns, Int. J. for Multiscale 
    Computational Engineering, 7(2)115-133(2009)
 
    Parameters set in such a way that the Figure 18 gets reproduced.
    At the moment no yarn damage assumed.
    '''
    from ibvpy.api import \
        TStepper as TS, RTraceGraph, RTraceDomainListField, TLoop, \
        TLine, IBVPSolve as IS, DOTSEval, BCSlice, FEGrid
    from ibvpy.fets.fets1D5.fets1D52l4uLRH import \
        FETS1D52L4ULRH,  MATS1DElastic, MATS1DPlastic, MATS1D5Bond
    from ibvpy.fets.fets1D5.fets1D52l6uLRH import FETS1D52L6ULRH
    from ibvpy.fets.fets1D5.fets1D52l8uLRH import FETS1D52L8ULRH
    from mathkit.mfn.mfn_line.mfn_line import \
        MFnLineArray
        
    from math import sqrt, pi as Pi

    # Concrete parameters
    A_concrete = 30. * 30.             # [mm^2] cross-sectional area of concrete
    E_concrete = 32000.0               # [MPa]  E-Modulus of concrete
    stiffness_concrete = E_concrete * A_concrete

    # Yarn parameters
    A_epoxy  = 1.760                   # [mm^2] cross-sectional area of epoxy
    A_glass  = 0.896                   # [mm^2] cross-sectional area of glass
    A_yarn   = A_epoxy + A_glass       # [mm^2] cross-sectional area of yarn
    P_yarn   = sqrt( 4 * A_yarn * Pi ) # [mm] perimeter of the yarn
    E_yarn   = 17000.0                 # [MPa]  effective E-Modulus of the impregnated yarn
    stiffness_fiber    = E_yarn * A_yarn
    
    # Bond parameters
    tau_max   = 13.0             # [N/mm^2] - frictional shear stress
    T_max     = tau_max * P_yarn # [N/mm] - frictional shear flow
    s_crit    = 0.028            # [mm] - onset of inelastic slip
    G         = T_max / s_crit   # [N/mm^2] - shear flow stiffness

    # Geometry
    L_e = 30.                            # [mm] - embedded length
    
    # Loading conditions
    u_max   = 0.1
    # f_max   = 1200
    
    # Material model construction
    mats_eval = MATS1D5Bond( mats_phase1 = MATS1DElastic( E = stiffness_fiber ),
                             mats_phase2 = MATS1DElastic( E = stiffness_concrete ),
                             mats_ifslip = MATS1DPlastic(E = G,
                                                         sigma_y = T_max,
                                                         K_bar = 0.,
                                                         H_bar = 0. ),
                             mats_ifopen = MATS1DElastic( E = 0. ))

    # Finite element construction
 #   fets_eval = FETS1D52L4ULRH( mats_eval = mats_eval ) # bilinear
 #   fets_eval = FETS1D52L6ULRH( mats_eval = mats_eval ) #quadratic
    fets_eval = FETS1D52L6ULRH( mats_eval = mats_eval ) #cubic
    # Dicretization
    domain = FEGrid( coord_max = (L_e, L_e/5.),
                     #shape   = (16,1), # for bilinear
                     #shape   = (8,1), # for quadratic
                     shape   = (4,1), # for cubic
                     fets_eval = fets_eval )

    end_dof = domain[-1,0,-1,0 ].dofs[0,0,0]
    ts = TS( dof_resultants = True,
         sdomain = domain,
         bcond_list =  [
                        # Matrix is fixed along the whole embedded length
                        BCSlice( var='u', value = 0., dims = [0], slice = domain[0,0,0,-1] ),
                                 
                        # Fixed fiber at the left end
                        BCSlice( var='u', value = 0., dims = [0],
                                  slice = domain[0,0,0,0] ),

                        # Fixed y-displacement - no crack opening
                        BCSlice( var='u', value = 0., dims = [1], slice = domain[:,:,:,:] ),
                                 
                        # Loading at the right end of the fiber
                        BCSlice( var='u', value = u_max, dims = [0], slice = domain[-1,0,-1,0] ) 
                        ],
         rtrace_list = [ RTraceGraph(name = 'Fi,right over u_right (iteration)' ,
                                     var_y = 'F_int', idx_y = end_dof,
                                     var_x = 'U_k', idx_x = end_dof),
                         RTraceDomainListField(name = 'slip', var = 'slip' ),
                         RTraceDomainListField(name = 'eps1', var = 'eps1' ),
                         RTraceDomainListField(name = 'eps2', var = 'eps2' ),
                         RTraceDomainListField(name = 'shear_flow', var = 'shear_flow' ),
                         RTraceDomainListField(name = 'sig1', var = 'sig1' ),
                         RTraceDomainListField(name = 'sig2', var = 'sig2' ),
                         RTraceDomainListField(name = 'Displacement', var = 'u' )                      
                ] )
    
    # Add the time-loop control
    tloop = TLoop( tstepper = ts, KMAX = 30, debug = False,
                   tline  = TLine( min = 0.0,  step = 0.1, max = 1.0 ))
    
    print tloop.eval()    
    # Put the whole stuff into the simulation-framework to map the
    # individual pieces of definition into the user interface.
    from ibvpy.plugins.ibvpy_app import IBVPyApp
    app = IBVPyApp( ibv_resource = tloop )
    app.main()

if __name__ == '__main__':
    run()