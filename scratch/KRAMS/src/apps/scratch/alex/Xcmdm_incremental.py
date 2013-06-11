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
# Created on May 29, 2009 by: rch


from enthought.traits.api import Float, Instance, Array, Int, Property, cached_property
from enthought.traits.ui.api import View, Item
from ibvpy.mats.mats2D.mats2D_explore import MATS2DExplore
from numpy import copy, array, hstack, loadtxt
from mathkit.mfn.mfn_line.mfn_line import MFnLineArray        
from mathkit.mfn.mfn_line.mfn_matplotlib_editor import MFnMatplotlibEditor
from mathkit.mfn.mfn_line.mfn_plot_adapter import MFnPlotAdapter
from scipy.optimize import brentq, newton, fsolve, brenth
from os.path import join

from ibvpy.core.scontext import SContext

# ---------------------------------------------------
# fitter: 
# ---------------------------------------------------

class MATS2DDamageFitter( MATS2DExplore ):
    '''
    Fitting algorithm for the damage function of the 
    quasi-ductile anisotropic material model.
    
    The algorithm uses the TLoop instance to proceed step 
    by step with the computation. The value of the damage function 
    for the time step t_n is identified iteratively by adjusting
    the values and evaluating the corresponding equilibrated stresses.
    
    The control parameters of the algorithm are:
    
    @param step_size: time step for fitting the damage parameter.
    @param tmax: end time for fitting, it might be also be set implicitly 
    for integrity = 1 - full damage of the material.
    '''

    # @todo: check if this yields different results 
    #        - check for pathologic behavior of the fitter:
    # is there a dependency on tmax or on n_steps??? 
    tmax      = 1.96666667e-03 
    
#    # @todo: after this value the linear elastic curve can't be fitted anymore
#    # because eps is greater then e_max_value
#    tmax      =     0.00167167 
    n_steps   = 10 
    step_size = tmax / n_steps 


    def run_through(self):
        '''Run the computation without fitting from the start to the end
        '''
        self.tloop.tline.max = self.tmax
        self.tloop.tline.step = self.step_size
        self.tloop.eval()
        print 'ending time',self.tloop.t_n1
        # show the response


    def run_step_by_step(self):
        '''Run the computation step by step from the start to the end
        '''
        n_steps = self.tmax / self.step_size
        self.tloop.tline.step = self.step_size
        current_time = 0.
        tmax = 0.
        
        for i in range( n_steps ):
            print 'STEP', i
            self.run_one_step()


    def run_trial_step( self ):
        '''Run the computation one step starting from the
        current time t_n to iterate the value for phi_new 
        which gives a fit with macroscopic stress curve.
        NOTE: The trial step does not update 'U_n' or 't_n'!
        '''
        print 'in trial step: self.tloop.U_n', self.tloop.U_n
        if len(self.tloop.U_n) == 0:
            current_U_n = self.tloop.tstepper.new_cntl_var()
            print 'U_n = None: tloop.tstepper.new_cntl_var()', self.tloop.tstepper.new_cntl_var()
        else:
            current_U_n = self.tloop.U_n[:]
        current_time = self.tloop.t_n
        self.run_one_step()

        # reset the current time back
        self.tloop.t_n = current_time
        self.tloop.U_n[:] = current_U_n[:]
        print '--- in trial step --- ' 

        self.tloop.tstepper.sctx.update_state_on = False


    def run_one_step(self):
        '''Run the computation one step starting from the
        current time t_n with the iterated value for phi_new
        in order to update TLoop and save the new phi value in
        the array ydata of PhiFnGeneral
        NOTE: The calculated step does update 'U_n' or 't_n'!
        '''
        self.tloop.tline.step = self.step_size
        current_time = self.tloop.t_n
        tmax = current_time + self.step_size
        self.tloop.tline.max = tmax
        self.tloop.eval()


    def get_target_data_tensile_test(self):
        '''subsidary method to read in the experimental 
        data (=target stress-strain curve for fitting) 
        and return it as numpy arrays
        '''
        # strains [-]:
        strain_arr = loadtxt( join( 'target_data', 'experimental_data', 'tensile_test_strain.dat') ) / 1000.
        # force [kN]: 
        force_arr = loadtxt( join( 'target_data', 'experimental_data', 'tensile_test_force.dat') )
        # cross-sectional area [mm**2]:
        A = 40.*100. 
        # stresses [MPa]: 
        stress_arr = force_arr / A * 1000.
        xdata_target = strain_arr[:]
        ydata_target = stress_arr[:]
        # make sure that the output arrays have the same number of values:
        data_len = min( len( xdata_target ), len( ydata_target ) ) 
        return xdata_target[:data_len], ydata_target[:data_len]
    
    
    #--------------------------------------------------
    # fitting data (from the macro-model (for verification))
    #--------------------------------------------------
    # specify target data:
#    target_data_spec = 'tensile_test'
    target_data_spec = 'linear_elastic'
#    target_data_spec = 'quasi_brittle'
#    target_data_spec = 'tension_stiffening'

    target_data = Property()
    @cached_property
    def _get_target_data(self):
        if self.target_data_spec == 'linear_elastic':
            from target_data.numerical_data.target_data_linear_elastic import xdata, ydata
        elif self.target_data_spec == 'quasi_brittle':
            from target_data.numerical_data.target_data_quasi_brittle import xdata, ydata
        elif self.target_data_spec == 'tension_stiffening':
            from target_data.numerical_data.target_data_tension_stiffening import xdata, ydata
        elif self.target_data_spec == 'tensile_test':
            xdata, ydata = self.get_target_data_tensile_test()
        print 'Selected target data specified as: ', self.target_data_spec
        print 'Selected target data: ydata:', ydata
        return xdata, ydata 
    
    xdata_target = Array
    def _xdata_target_default(self):
        return self.target_data[0]
    
    ydata_target = Array
    def _ydata_target_default(self):
        return self.target_data[1]

    #--------------------------------------------------
    # interpolation function for fitting data:
    #--------------------------------------------------
    mfn_line_array_target = Instance( MFnLineArray )
    def _mfn_line_array_target_default(self):
        mfn = MFnLineArray( xdata = self.xdata_target, ydata = self.ydata_target )
        mfn.data_changed = True
        return mfn
    
    
    def init(self):
        #--------------------------------------------------
        # for fitting use 'General'-function for 'phi_fn': 
        #--------------------------------------------------
        # The value pair for the piecewise linear definition
        # of 'phi_fn' value consists of current strain and the
        # iterated 'phi_value'. The microplanes with a lower
        # microplane strain level use an interpolated value 
        # for 'phi' 
        self.fitted_phi_fn = self.mats2D_eval.polar_fn.phi_fn.mfn
        self.fitted_phi_fn.xdata = [0]
        self.fitted_phi_fn.ydata = [1]
        self.fitted_phi_fn.data_changed = True
#        self.mats2D_eval.polar_fn.phi_fn.refresh_plot()
#        print 'value', self.mats2D_eval.polar_fn.phi_fn.mfn.get_value(0.5)
            
            
    def get_lack_of_fit( self, phi_trial ):
        '''Return the difference between the macroscopic stress calculated
        based on the value of phi_trial (damage at the next step) and the
        macroscopic stress defined as target data (=fitting curve)
        '''
        # current maximum macroscopic strain corresponds to control variable
        current_time = self.tloop.t_n
        print '    current_time = ', current_time    
        print '    step_size = ', self.step_size    

        # get the maximum microplane strain:
        if len(self.tloop.U_n) == 0:
            e_max_value = 0.0
        else:
            U_n = copy( self.tloop.U_n[:] )
            print '    U_n = ', U_n    

            d_U = copy( self.tloop.d_U[:] )
            print '    d_U = ', d_U    
            
            U_n1 = U_n + d_U 
            print '    U_n1 = ', U_n1    
            
            sctx = self.tloop.tstepper.sctx
            # do not update the state array in the trial steps
            sctx.update_state_on = False
            e_max_arr = self.mats2D_eval._get_state_variables( sctx, U_n1 )
            e_max_value = max( e_max_arr )
            print 'e_max_arr = ', e_max_arr    

            
        print '    e_max_value = ', e_max_value    
#        print 'current_time = ', current_time    
                
        # use old phi_value within trial step:
        x = hstack([ self.fitted_phi_fn.xdata[:], current_time + self.step_size ])
#        x = hstack([ self.fitted_phi_fn.xdata[:], e_max_value ])
        y = hstack([ self.fitted_phi_fn.ydata[:], phi_trial     ])
        print 'in get_lack_of_fit: current time', current_time + self.step_size
        self.fitted_phi_fn.set( xdata = x, ydata = y )
        self.fitted_phi_fn.data_changed = True
                    
        # try the next equilibrium
        self.run_trial_step()

        # remove the trial value from the damage function:
        x = self.fitted_phi_fn.xdata[:-1]
        y = self.fitted_phi_fn.ydata[:-1]
        self.fitted_phi_fn.set( xdata = x, ydata = y )
        self.fitted_phi_fn.data_changed = True
        
        # compare the obtained stress with the measured response
            
        # get calculated value for sig_app based
        # on the current value of 'phi_trial':
        self.tloop.rtrace_mngr.rtrace_bound_list[1].redraw()
        sig_app_trial = self.tloop.rtrace_mngr.rtrace_bound_list[1].trace.ydata[-1]
#        print 'whole sig_app', self.tloop.rtrace_mngr.rtrace_bound_list[1].trace.ydata
        # get corresponding value from the fitted data:
        sig_app_target = self.mfn_line_array_target.get_value( current_time + self.step_size )
        lack_of_fit = sig_app_trial - sig_app_target

        print '    phi_trial:    ', phi_trial 
        print '    sig_app_trial ', sig_app_trial    
        print '    sig_app_target', sig_app_target

        return lack_of_fit
        
        
    def fit_response(self):
        '''iterate phi_trian in each incremental step such that the
        lack of fit between the calculated stress and the fitting curve vanishes
        '''
        phi_old = 1.0
        
        # @todo: is there a possibility to memorize the trials 
        # of 'brentq' and use the last best value instead of phi_old?
        phi_trial_list = []
        
        for n in range( self.n_steps ):

            # use scipy-functionality to get the iterated value of phi_new
            # If both phi_trial = phi_old and phi_trial = 0. return the same 
            # lack of fit no sign change is given which is a requirement
            # for the function call 'brentq'. In this case use old value
            # for phi_trial is used and tloop moves on one step 

            phi_new = phi_old
            try:
                # @todo: check if there should be optional arguments for method
                # 'brentq' such as 'xtol' or 'maxiter' (so far default values ar used)
                phi_new = brentq( self.get_lack_of_fit, 0., phi_old )#, xtol = 1e-3)

                # @todo: check if this is gives better fitting results; faster? 
#                phi_new = brenth( self.get_lack_of_fit, 0., phi_old )
            except ValueError:
                 print "No sign change between get_lack_of_fit(phi_old) and get_lack_of_fit(0.). \n Use old value for phi_trial. phi_old = ", phi_old
                             
            # current time corresponds to the current strain applied
            current_time = self.tloop.t_n
    
            # replace old 'phi_value' with iterated value:
            phi_old = phi_new
                        
         
#            
#            # get the maximum microplane strain:
#            if len(self.tloop.U_n) == 0:
#                e_max_value = 0.0
#            else:
#                U_n = copy( self.tloop.U_n[:] )
#                print 'U_n = ', U_n    
#    
#                d_U = copy( self.tloop.d_U[:] )
#                print 'd_U = ', d_U    
#                
#                U_n1 = U_n + d_U 
#                print 'U_n1 = ', U_n1    
#                
#                sctx = self.tloop.tstepper.sctx
#                # do not update the state array in the trial steps
#                sctx.update_state_on = False
#                e_max_arr = self.mats2D_eval._get_state_variables( sctx, U_n1 )
#                e_max_value = max( e_max_arr )
#                print 'e_max_arr = ', e_max_arr                            
                        
                        
            # update phi_data:
            x = hstack([ self.fitted_phi_fn.xdata[:], current_time + self.step_size  ])
#            x = hstack([ self.fitted_phi_fn.xdata[:], e_max_value  ])
            y = hstack([ self.fitted_phi_fn.ydata[:],  phi_new    ])
            self.fitted_phi_fn.set( xdata = x, ydata = y )
            self.fitted_phi_fn.data_changed = True
            
            # @todo: (remove after testing):
            # print out the appended data for the phi function 
            # but only for the first 4 converged steps:
            if len(x)<=4:
                print 'update xdata: in step', n+1, ': ', x
                print 'update ydata: in step', n+1, ': ', y
            
            # run one step with the iterated value for phi in order to
            # update the state array and to move forward one step:
            print '### run_one_step ###: current time: ', self.tloop.t_n
            self.run_one_step()


     # used to plot the target curve:
#    target_view = View( Item( 'mfn_line_array_target', show_label = False,
#                              editor = MFnMatplotlibEditor( adapter = MFnPlotAdapter( var_x = 'xdata_target',
#                                                                                      var_y = 'ydata_target' ))),
#                        height = 0.7, width = 0.8, resizable = True, buttons = ['Ok', 'Cancel'] )
                        


if __name__ == '__main__':
    #--------------------------------------------------------------------------------
    # Example using the mats2d_explore 
    #--------------------------------------------------------------------------------
    from ibvpy.mats.mats2D.mats2D_explore import MATS2DExplore
    from ibvpy.mats.mats2D.mats2D_rtrace_cylinder import MATS2DRTraceCylinder
    from ibvpy.mats.mats2D.mats2D_cmdm.mats2D_cmdm_rtrace_Gf_mic import MATS2DMicroplaneDamageTraceGfmic,\
        MATS2DMicroplaneDamageTraceEtmic, MATS2DMicroplaneDamageTraceUtmic
    from ibvpy.mats.mats2D.mats2D_cmdm.mats2D_cmdm_rtrace_Gf_mac import MATS2DMicroplaneDamageTraceGfmac,\
        MATS2DMicroplaneDamageTraceEtmac, MATS2DMicroplaneDamageTraceUtmac
    from ibvpy.mats.mats2D.mats2D_cmdm.mats2D_cmdm import MA2DMicroplaneDamage   
    from ibvpy.api import RTraceGraph, RTraceArraySnapshot 
    from ibvpy.mats.mats2D.mats2D_cmdm.mats2D_cmdm_damage_fn import PhiFnQuasiBrittle
    
    from mathkit.mfn.mfn_line.mfn_line import MFnLineArray        
    from numpy import array, hstack    
    

    fitter = MATS2DDamageFitter \
            (          mats2D_eval = MA2DMicroplaneDamage( elastic_debug = False,
                                                           polar_fn_class = 'Isotropic polar function'  
                        ),
                       rtrace_list = [ 
                                       MATS2DRTraceCylinder(name = 'Laterne',
                                                            var_axis    = 'time', idx_axis = 0,
                                                            var_surface = 'microplane_damage',
                                                            update_on = 'update' ),
                                                            
                                       RTraceGraph(name = 'strain - stress',
                                                   var_x = 'eps_app', idx_x = 0,
                                                   var_y = 'sig_app', idx_y = 0,
                                                   record_on = 'update' ),
                                       RTraceGraph(name = 'strain - strain',
                                                   var_x = 'eps_app', idx_x = 0,
                                                   var_y = 'eps_app', idx_y = 1,
                                                   record_on = 'update' ),
                                       RTraceGraph(name = 'stress - stress',
                                                   var_x = 'sig_app', idx_x = 0,
                                                   var_y = 'sig_app', idx_y = 1,
                                                   record_on = 'update' ),

#                                       RTraceGraph(name = 'strain - stress',
#                                                   var_x = 'eps_app_v3d', idx_x = 0,
#                                                   var_y = 'sig_app_v3d', idx_y = 0,
#                                                   record_on = 'update' ),
#                                       RTraceGraph(name = 'strain - strain',
#                                                   var_x = 'eps_app_v3d', idx_x = 0,
#                                                   var_y = 'eps_app_v3d', idx_y = 1,
#                                                   record_on = 'update' ),
#                                       RTraceGraph(name = 'stress - stress',
#                                                   var_x = 'sig_app_v3d', idx_x = 0,
#                                                   var_y = 'sig_app_v3d', idx_y = 1,
#                                                   record_on = 'update' ),

                                       RTraceGraph(name = 'time - microplane damage',
                                                   var_x = 'time', idx_x = 0,
                                                   var_y = 'microplane_damage', idx_y = 0,
                                                   update_on = 'update' ),

                                       # e_equiv, s_equiv
                                       RTraceGraph(name = 'e_equiv - s_equiv',
                                                   var_x = 'e_equiv_arr', idx_x = 0,
                                                   var_y = 's_equiv_arr', idx_y = 0,
                                                   update_on = 'update' ),
                                    
                                       # e_N, s_N:            
                                       RTraceGraph(name = 'e_N - s_N',
                                                   var_x = 'e_N_arr', idx_x = 0,
                                                   var_y = 's_N_arr', idx_y = 0,
                                                   update_on = 'update' ),

                                       # e_T, s_T:            
                                       RTraceGraph(name = 'e_T - s_T',
                                                   var_x = 'e_T_arr', idx_x = 0,
                                                   var_y = 's_T_arr', idx_y = 0,
                                                   update_on = 'update' ),

                                       RTraceArraySnapshot(name = 'equiv_projection',
                                                           var = 'equiv_projection',
                                                           record_on = 'update' ),
                                                   
                                       RTraceArraySnapshot(name = 'microplane damage',
                                                           var = 'microplane_damage',
                                                           record_on = 'update' ),

                                       RTraceArraySnapshot(name = 'e_equiv',
                                                           var = 'e_equiv_arr',
                                                           record_on = 'update' ),
                                       RTraceArraySnapshot(name = 's_equiv',
                                                           var = 's_equiv_arr',
                                                           record_on = 'update' ),
                                       RTraceArraySnapshot(name = 'e_N',
                                                           var = 'e_N_arr',
                                                           record_on = 'update' ),
                                       RTraceArraySnapshot(name = 's_N',
                                                           var = 's_N_arr',
                                                           record_on = 'update' ),
                                       RTraceArraySnapshot(name = 'e_T',
                                                           var = 'e_T_arr',
                                                           record_on = 'update' ),
                                       RTraceArraySnapshot(name = 's_T',
                                                           var = 's_T_arr',
                                                           record_on = 'update' ),
                                                           
                                     ]
                       )
            
    #---------------------------
    # Default settings:
    #---------------------------

    fitter.mats2D_eval.polar_fn.n_mp = 7
    print 'n_mp:          ', fitter.mats2D_eval.polar_fn.n_mp

    fitter.alpha_degree = 0.
    print 'alpha_degree:  ', fitter.alpha_degree

    fitter.mats2D_eval.stress_state = 'plane_stress'
#    fitter.mats2D_eval.stress_state = 'plane_strain'

    fitter.mats2D_eval.symmetrization = 'sum-type'
    print 'symmetrization:', fitter.mats2D_eval.symmetrization

    fitter.mats2D_eval.model_version = 'compliance'
#    fitter.mats2D_eval.model_version = 'stiffness'
    print 'model_version: ', fitter.mats2D_eval.model_version

    # Default settings for 'polar_fn_class':
    fitter.mats2D_eval.polar_fn_class = 'Isotropic polar function'

    # Default settings of 'PolarFnBase' for 'phi_fn_class':
    fitter.mats2D_eval.polar_fn.phi_fn_class = 'General'

    #---------------------------
    # run fitter:
    #---------------------------
    fitter.init()
    fit_response = fitter.fit_response()





    #---------------------------
    # basic testing of fitter methods:
    #---------------------------

    basic_tests = False
    if basic_tests:        
        fitter.run_through()
        #    fitter.tloop.rtrace_mngr.rtrace_bound_list[1].configure_traits()
        fitter.tloop.rtrace_mngr.rtrace_bound_list[1].redraw()
        last_strain_run_through = fitter.tloop.rtrace_mngr.rtrace_bound_list[1].trace.xdata[:]
        last_stress_run_through = fitter.tloop.rtrace_mngr.rtrace_bound_list[1].trace.ydata[:]
        print 'last strain (run-through) value' , last_strain_run_through
        print 'last stress (run-through) value' , last_stress_run_through
    
        fitter.tloop.reset()
        fitter.run_step_by_step()
        #fitter.tloop.rtrace_mngr.rtrace_bound_list[1].configure_traits()
        fitter.tloop.rtrace_mngr.rtrace_bound_list[1].redraw()
        last_strain_step_by_step = fitter.tloop.rtrace_mngr.rtrace_bound_list[1].trace.xdata[:]
        last_stress_step_by_step = fitter.tloop.rtrace_mngr.rtrace_bound_list[1].trace.ydata[:]
        print 'last stress (step-by-step) value', last_stress_step_by_step
        
        fitter.run_trial_step()
        fitter.run_trial_step()
        fitter.tloop.rtrace_mngr.rtrace_bound_list[1].redraw()
        strain_after_trial_steps = fitter.tloop.rtrace_mngr.rtrace_bound_list[1].trace.xdata[:]
        stress_after_trial_steps = fitter.tloop.rtrace_mngr.rtrace_bound_list[1].trace.ydata[:]
        print 'stress after trial', stress_after_trial_steps
    
        fitter.init()
        #fitter.mats2D_eval.configure_traits()
        lof = fitter.get_lack_of_fit( 1.0 )
        print '1',lof
        lof = fitter.get_lack_of_fit( 0.9 )
        print '2',lof
    
        #fitter.tloop.rtrace_mngr.configure_traits()
        fitter.run_trial_step()
    
    
    
    from ibvpy.plugins.ibvpy_app import IBVPyApp
    ibvpy_app = IBVPyApp( ibv_resource = fitter )
    ibvpy_app.main()       
