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
# Created on Jul 29, 2009 by: rchx

#----------------------- example --------------------


if __name__ == '__main__':

    from math import fabs
    from numpy import loadtxt
    from os.path import join
        
    from ibvpy.fets.fets2D5.fets2D58h import FETS2D58H
    from ibvpy.fets.fets2D5.fets2D58h20u import FETS2D58H20U
    from ibvpy.bcond.bc_slice import BCSlice    
    from ibvpy.api import \
        TStepper as TS, RTraceGraph, RTraceDomainListField, TLoop, \
        TLine, BCDofGroup, IBVPSolve as IS, DOTSEval
    from ibvpy.fets.fets3D.fets3D8h import FETS3D8H
    from ibvpy.fets.fets3D.fets3D8h20u import FETS3D8H20U    
    from ibvpy.mats.mats3D.mats3D_elastic.mats3D_elastic import MATS3DElastic
    from ibvpy.mats.mats2D.mats2D_elastic.mats2D_elastic import MATS2DElastic

    from mathkit.mfn.mfn_line.mfn_line import MFnLineArray

    from ibvpy.mats.mats2D5.mats2D5_cmdm.mats2D5_cmdm import \
        MATS2D5MicroplaneDamage
        
    from ibvpy.mats.matsXD.matsXD_cmdm.matsXD_cmdm_phi_fn import \
        PhiFnStrainSoftening, PhiFnGeneral, PhiFnStrainHardening, \
        PhiFnStrainHardeningLinear, PhiFnGeneralExtended    

    from ibvpy.mats.matsXD.matsXD_cmdm.matsXD_cmdm_phi_fn import \
        PhiFnGeneralExtended
            
    from promod.simdb import \
        SimDB
    
    simdb = SimDB()
    
    import pickle
    
            
 
    from promod.matdb.trc.ccs_unit_cell import CCSUnitCell
    
    #ccsuc = CCSUnitCell.db['PZ-0708-1_MAG-07-03_000300_90_0']
    ccsuc = CCSUnitCell.db['FIL-10-09_2D-02-06a_0.00273_90_0']
    mfn = ccsuc.damage_function_list[0].damage_function
    
    phi_fn = PhiFnGeneralExtended( mfn = mfn )
    
    rho = 0.0321 # corresponds to 'TT11-10a_2D-02-06a.mats' (E_c =  
    E_f = 70000
    E_m = 28700  
    E_c = E_m * ( 1 - rho ) + E_f * rho


    #phi_fn = PhiFnStrainHardening(Epp = 0.00025, Efp =0.0003, Dfp = 0.35, Elimit = 0.008 )          
    #phi_fn = PhiFnStrainHardening(Epp = 0.0001, Efp =10, Dfp = 0.35, Elimit = 0.1 )
#        phi_fn = PhiFnStrainHardening(Epp = 0.0001, Efp =0.01, Dfp = 0.4, Elimit = 0.006 )
#
#        Epp = 0.0001
#        Dfp = 0.4
#        Efr = 0.0
#        Efp = - ( -1 + Dfp ) * Epp / ( Dfp - Efr ) 
#        phi_fn = PhiFnStrainHardening(Epp = Epp, Efp = Efp, Dfp = Dfp, Elimit = 0.006 )        

#        phi_fn = PhiFnStrainHardeningLinear( E_m = E_m, E_f = E_f, rho = rho,
#                                             sigma_0 = 5.0,  
#                                             alpha = 0.5, beta = 0.5, Elimit = 0.06 )
    
    print 'C_c', E_c
    mats_eval = MATS2D5MicroplaneDamage( E = E_c, nu = 0.2,
                                         n_mp = 10,
                                         symmetrization = 'sum-type',
                                         model_version = 'compliance',
                                         phi_fn = phi_fn ) 

    fets_eval = FETS2D58H20U(mats_eval  = mats_eval)

    #fets_eval.vtk_r *= 0.577
                                        
    from ibvpy.mesh.fe_grid import FEGrid
    
    L = 1.25 # [m]
    L_quater = L / 2. # [m]
    thickness = 0.03

    w_max = - 0.05 # [m]
    
    n_xy = 8
    n_z  = 2

    # Discretization
    domain = FEGrid( coord_max = (L_quater,L_quater,thickness ), 
                     shape   = (n_xy, n_xy, n_z),
                     fets_eval = fets_eval)

    upper_left_corner = domain[-1,-1,-1,-1,-1,-1].dofs[0,0,2]

    bc_symplane_yz  = BCSlice( var = 'u', value = 0., dims = [0], slice = domain[-1,:,:,-1,:,:] )
    bc_symplane_xz  = BCSlice( var = 'u', value = 0., dims = [1], slice = domain[:,-1,:,:,-1,:] )
    bc_support_000  = BCSlice( var = 'u', value = 0., dims = [2], slice = domain[0,0,0,0,0,0] )
    bc_center_w     = BCSlice( var = 'u', value = w_max, dims = [2], slice = domain[-1,-1,-1,-1,-1,-1] )

    f_w_diagram = RTraceGraph(name = 'time - reaction 2',
                                       var_x = 'U_k', idx_x = upper_left_corner,
                                       var_y = 'F_int', idx_y = 2,
                                       record_on = 'update',
                                       transform_x = '-x * 1000', # %g * x' % ( fabs( w_max ),),
                                       transform_y = '4 * 1000 * y' )                                     
    ts = TS(
            sdomain = domain,
            bcond_list = [bc_symplane_yz, bc_symplane_xz, bc_support_000, 
                          bc_center_w,
                        ],
             rtrace_list = [      
                         f_w_diagram,
                         RTraceDomainListField(name = 'Stress' ,
                                        var = 'sig_app', idx = 0, warp = True, 
                                        record_on = 'update'),
                         RTraceDomainListField(name = 'Strain' ,
                                        var = 'eps_app', idx = 0, warp = True, 
                                        record_on = 'update'),                        
                         RTraceDomainListField(name = 'Damage' ,
                                        var = 'omega_mtx', idx = 0, warp = True, 
                                        record_on = 'update'),
                         RTraceDomainListField(name = 'IDamage' ,
                                        position = 'int_pnts',
                                        var = 'omega_mtx', idx = 0, 
                                        record_on = 'update'),                                        
                                        ]             
            )
    
    # Add the time-loop control

    #            tloop = TLoop( tstepper = ts, 
    #                           KMAX = 50, 
    #                           RESETMAX = 0, 
    #                           tolerance = 5e-4, 
    #                           debug = False,
    #                           tline  = TLine( min = 0.0,  step = 0.1, max = 1.0 ))
    
    tloop = TLoop( tstepper = ts,
           KMAX = 20, 
           RESETMAX = 0, 
           tolerance = 0.01, # 5e-4, 
           debug = False,
           tline  = TLine( min = 0.0,  step = 0.02, max = 1.0 ) )                   

    from pickle import dump, load

    do = 'show_last_result'
#    do = 'eval'
    
    if do == 'eval':
        print 'dof_number', upper_left_corner
        U = tloop.eval()
        #print 'U', U[ upper_left_corner ]
    
        f_w_diagram.refresh()
        file = open( 'f_w_diagram.pickle', 'w' )
        dump( f_w_diagram.trace, file )
        file.close()

    if do == 'show_last_result':
        from promod.exdb.ex_run import ExRun
        import pylab as p

#        f_w_diagram.trace.mpl_plot( p, color = 'red' )

        file = open( 'f_w_diagram.pickle', 'r' )
        trace = load( file )
        p.plot( trace.xdata, trace.ydata, color = 'blue' )

        file = open( 'f_w_diagram_TT11-10a.pickle', 'r' )
        trace = load( file )
        p.plot( trace.xdata, trace.ydata, color = 'red' )
        
#        file = open( 'f_w_diagram.pickle_SimModel', 'r' )
#        trace = load( file )
#        p.plot( trace.xdata, trace.ydata, color = 'black' )
#        
        file = open( 'f_w_diagram_SimModel.pickle', 'r' )
        trace = load( file )
        p.plot( trace.xdata, trace.ydata, color = 'green' )

#        path = join( simdb.exdata_dir, 'plate_tests', 'PT-9a' )
#        tests = [ 'PT01-9a.DAT', 'PT02-9a.DAT' , 'PT09-9a.DAT' ]

        path = join( simdb.exdata_dir, 'plate_tests', 'PT-10a' )
        tests = [ 'PT10-10a.DAT', 'PT11-10a.DAT' , 'PT12-10a.DAT' ]
        
        for t in tests:
            ex_path = join( path, t )        
            ex_run = ExRun( ex_path )

            ex_run.ex_type._plot_smoothed_force_deflection_center( p )
                    
        p.show()
    
    # Put the whole stuff into the simulation-framework to map the
    # individual pieces of definition into the user interface.
    if do == 'ui':
        from ibvpy.plugins.ibvpy_app import IBVPyApp
        app = IBVPyApp( ibv_resource = tloop )
        app.main()
