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


from enthought.traits.api import \
    Float, Instance, Array, Int, Property, cached_property, on_trait_change, Bool, \
    HasTraits
from enthought.traits.ui.api import View, Item
from ibvpy.mats.mats_explore import MATSExplore
from ibvpy.mats.mats2D.mats2D_explore import MATS2DExplore
from numpy import copy, array, hstack, loadtxt, savetxt
from mathkit.mfn.mfn_line.mfn_line import MFnLineArray        
from mathkit.mfn.mfn_line.mfn_matplotlib_editor import MFnMatplotlibEditor
from mathkit.mfn.mfn_line.mfn_chaco_editor import MFnChacoEditor

from mathkit.mfn.mfn_line.mfn_plot_adapter import MFnPlotAdapter
from scipy.optimize import brentq, newton, fsolve, brenth
from os.path import join
from ibvpy.core.tloop import TLoop, TLine
#from tloop import TLoop, TLine
from ibvpy.core.scontext import SContext
from ibvpy.core.tstepper import TStepper

from promod.exdb.ex_run import ExRun

# ---------------------------------------------------
# fitter: 
# ---------------------------------------------------

class MATSDamageFitter( MATSExplore ):
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

    #max_load  = 0.2*9.7599/ 1000.  #1.96666667e-03 

    max_load = Property( Float )
    def _get_max_load(self):
        return self.xdata_target[-1] * 1.1
                         
    n_steps   = 30
    
    step_size = Property( Float )
    def _get_step_size(self):
        return self.max_load / self.n_steps
    
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
        print '--------- run trial step: --------- '
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
        print '--------- end of trial step --------- '

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
        self.update_e_max_value_new = True

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
    
    def get_target_data_tsp(self):
        '''subsidary method to read in the experimental 
        data (=target stress-strain curve for fitting) 
        and return it as numpy arrays
        '''
        # strains [-]:
        strain_arr = loadtxt( join( 'target_data', 'experimental_data', 'tsp_03_eps.dat') ) / 1000.
        # stresses [MPa]: 
        stress_arr = loadtxt( join( 'target_data', 'experimental_data', 'tsp_03_sig.dat') )
        xdata_target = strain_arr[:]
        ydata_target = stress_arr[:]
        # make sure that the output arrays have the same number of values:
        data_len = min( len( xdata_target ), len( ydata_target ) ) 
        return xdata_target[:data_len], ydata_target[:data_len]
    
    ex_run = Instance( ExRun )
    def get_target_data_exdb_tensile_test(self):
        '''Use the data from the ExDB
        '''
        return self.ex_run.ex_type.eps_fit, self.ex_run.ex_type.sig_c_fit
    #--------------------------------------------------
    # fitting data (from the macro-model (for verification))
    #--------------------------------------------------
    # specify target data:
#    target_data_spec = 'tensile_test'
#    target_data_spec = 'linear_elastic'
#    target_data_spec = 'quasi_brittle'
#    target_data_spec = 'tension_stiffening'
    target_data_spec = 'ex_tensile_test'
#    target_data_spec = 'tsp'

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
        elif self.target_data_spec == 'tsp':
            xdata, ydata = self.get_target_data_tsp()
        elif self.target_data_spec == 'ex_tensile_test':
            xdata, ydata = self.get_target_data_exdb_tensile_test()
        print 'Selected target data specified as: ', self.target_data_spec
        return xdata, ydata 
    
    xdata_target = Property( Array )
    def _get_xdata_target(self):
        return self.target_data[0]
    
    ydata_target = Property( Array )
    def _get_ydata_target(self):
        return self.target_data[1]

    #--------------------------------------------------
    # interpolation function for fitting data:
    #--------------------------------------------------
    mfn_line_array_target = Instance( MFnLineArray )
    def _mfn_line_array_target_default(self):
        mfn = MFnLineArray( xdata = self.xdata_target, ydata = self.ydata_target )
        mfn.data_changed = True
        return mfn
    
    fitted_phi_fn = Instance( MFnLineArray )
    
    def init(self):
        #--------------------------------------------------
        # for fitting use 'General'-function for 'phi_fn': 
        #--------------------------------------------------
        # The value pair for the piecewise linear definition
        # of 'phi_fn' value consists of current strain and the
        # iterated 'phi_value'. The microplanes with a lower
        # microplane strain level use an interpolated value 
        # for 'phi' 
        self.fitted_phi_fn = self.dim.mats_eval.phi_fn.mfn
        self.fitted_phi_fn.xdata = [0]
        self.fitted_phi_fn.ydata = [1]
        self.fitted_phi_fn.data_changed = True
        # initialize TLoop parameters:
        self.tloop.setup()

    def get_lack_of_fit( self, phi_trial ):
        '''Return the difference between the macroscopic stress calculated
        based on the value of phi_trial (damage at the next step) and the
        macroscopic stress defined as target data (=fitting curve)
        '''
        print '\n'
        print "#'get_lack_of_fit' for the trial value "
        print '    phi_trial    = ', phi_trial 

        # value of the principle macroscopic strain corresponds to control variable
        current_time = self.tloop.t_n
        print '    current_time = ', current_time    
        print '    step_size    = ', self.step_size    

        # ------------------------------------                
        # add new pair in fitted_phi_fn 
        # ------------------------------------                
        # consisting of 'e_max_value_new' and 'phi_trial'
        x = hstack([ self.fitted_phi_fn.xdata[:], current_time + self.step_size ])
        y = hstack([ self.fitted_phi_fn.ydata[:], phi_trial ])
        self.fitted_phi_fn.set( xdata = x, ydata = y )
        self.fitted_phi_fn.data_changed = True
                    
        # ------------------------------------                
        # get state array before trial:
        # ------------------------------------                
        mats_state_array_old = copy( self.tloop.tstepper.sctx.mats_state_array )
                   
        # ------------------------------------                
        # run trial step: 
        # ------------------------------------                
        print '    reset current_U_n   =', self.tloop.U_n           
        print 'CURRENT PHI', self.dim.mats_eval.phi_fn.mfn.ydata             
        # try the next equilibrium
        self.run_trial_step()

        # ------------------------------------                
        # reset mats_state_array:
        # ------------------------------------                
        # Note: the material state array (i.e. the maximum microstrains) are 
        # updated within the iterations of each trial step, therefore a reset
        # is necessary in order to start each trial step with the same state variables  
        self.tloop.tstepper.sctx.mats_state_array[:] = mats_state_array_old[:]
        print '    reset state array'

        # ------------------------------------                
        # remove trial value in fitted_phi_fn 
        # ------------------------------------                
        x = self.fitted_phi_fn.xdata[:-1]
        y = self.fitted_phi_fn.ydata[:-1]
        self.fitted_phi_fn.set( xdata = x, ydata = y )
        self.fitted_phi_fn.data_changed = True
        
        # ------------------------------------                
        # get the lack of fit
        # ------------------------------------                
        # get calculated value for sig_app based on the current value of 'phi_trial':
        # and evaluate the difference between the obtained stress and the measured response
        self.tloop.rtrace_mngr.rtrace_bound_list[0].redraw()
        sig_app_trial = self.tloop.rtrace_mngr.rtrace_bound_list[0].trace.ydata[-1]
        # get corresponding value from the fitted data:
        sig_app_target = self.mfn_line_array_target.get_value( current_time + self.step_size )
        # absolut error:
        lack_of_fit_absolut = sig_app_trial - sig_app_target
        # relative error:
        lack_of_fit_relative = lack_of_fit_absolut / sig_app_target

        print '    sig_app_trial ', sig_app_trial    
        print '    sig_app_target', sig_app_target
        print '    lack_of_fit_absolut   ', lack_of_fit_absolut
        print '    lack_of_fit_relative  ', lack_of_fit_relative
        print '# get_lack_of_fit # END '
        return lack_of_fit_relative
        
        
    def fit_response(self):
        '''iterate phi_trial in each incremental step such that the
        relative lack of fit between the calculated stress and the target
        curve is smaller then xtol defined in function 'brentq'
        '''
        phi_old = 1.0
        
        for n in range( self.n_steps ):

            phi_new = phi_old
            
            # use scipy-functionality to get the iterated value of phi_new
            # If the trial value calculated with phi_trial = phi_old
            # is smaller then target value get_lack_of_fit has no sign change
            # for phi_trial = phi_old and phi_trial = 0. which is a requirement
            # for the function call 'brentq'. In this case the old value
            # for phi_trial is used and tloop moves on one step 
            try:
                # The method brentq has optional arguments such as
                #   'xtol'    - absolut error (default value = 1.0e-12)
                #   'rtol'    - relative error (not supported at the time)
                #   'maxiter' - maximum numbers of iterations used
                #
                # Here xtol is used to specify the allowed RELATIVE error!
                # therefore the relative lack of fit is returned in 
                # method 'get_lack_of_fit' 
                _xtol = 1.0e-6
                phi_new = brentq( self.get_lack_of_fit, 0., phi_old, xtol = _xtol)
                # @todo: check if this is gives better fitting results; faster? 
#                phi_new = brenth( self.get_lack_of_fit, 0., phi_old )
            except ValueError:              
                a = self.get_lack_of_fit( 0. )    
                b = self.get_lack_of_fit( phi_old )    
                print 'No sign change between get_lack_of_fit(phi_old) = ', a, ' and ' 
                print 'get_lack_of_fit(0.) = ', b
                print 'Use old value for phi_trial. phi_old = ', phi_old
                             
            # current time corresponds to the current strain applied
            current_time = self.tloop.t_n
    
            # replace old 'phi_value' with iterated value:
            phi_old = phi_new
                        
            # get mats_state_array:
            mats_state_array = copy( self.tloop.tstepper.sctx.mats_state_array )
    
            # update phi_data:
            x = hstack([ self.fitted_phi_fn.xdata[:], current_time + self.step_size  ])
            y = hstack([ self.fitted_phi_fn.ydata[:],          phi_new             ])
            self.fitted_phi_fn.set( xdata = x, ydata = y )
            self.fitted_phi_fn.data_changed = True

            # print out the appended data for the phi function 
            # (only for the first 4 converged steps):
            print 'len(x)', len(x)
            print 'n', n
            print 'self.n_steps', self.n_steps
            if n == (self.n_steps - 1):
                print 'update xdata: ', x
                print 'update ydata: ', y
                output_file = open('fitted_phi_fn.out','w')
                savetxt('fitted_phi_fn.out', (x,y))

            # run one step with the iterated value for phi in order to
            # update the state array and to move forward one step:
            print '\n'
            print '### run_one_step ###' 
            print '### step', n   ,'###'
            print '### current time:', current_time  
            self.run_one_step()
        self.fitted_phi_fn.changed = True

    # used to plot the target curve:
    traits_view = View( Item( 'mfn_line_array_target', show_label = False,
                              editor = MFnChacoEditor( adapter = MFnPlotAdapter())),
                        Item( 'fitted_phi_fn', show_label = False,
                              editor = MFnMatplotlibEditor( adapter = MFnPlotAdapter( var_x = 'xdata_target',
                                                                                      var_y = 'ydata_target' ))),
                        height = 0.7, width = 0.8, resizable = True, buttons = ['OK', 'Cancel'] )
                        
def run():
    #--------------------------------------------------------------------------------
    # Example using the mats2d_explore 
    #--------------------------------------------------------------------------------
    from ibvpy.mats.mats2D.mats2D_explore import MATS2DExplore
    from ibvpy.mats.mats2D.mats2D_rtrace_cylinder import MATS2DRTraceCylinder
    from ibvpy.mats.mats2D.mats2D_cmdm.mats2D_cmdm_rtrace_Gf_mic import \
        MATS2DMicroplaneDamageTraceGfmic,\
        MATS2DMicroplaneDamageTraceEtmic, MATS2DMicroplaneDamageTraceUtmic
    from ibvpy.mats.mats2D.mats2D_cmdm.mats2D_cmdm_rtrace_Gf_mac import \
        MATS2DMicroplaneDamageTraceGfmac,\
        MATS2DMicroplaneDamageTraceEtmac, MATS2DMicroplaneDamageTraceUtmac
    from ibvpy.mats.mats2D.mats2D_cmdm.mats2D_cmdm import \
        MATS2DMicroplaneDamage, PhiFnGeneral, PhiFnStrainSoftening   
    from ibvpy.api import RTraceGraph, RTraceArraySnapshot 
    
    from mathkit.mfn.mfn_line.mfn_line import MFnLineArray        
    from numpy import array, hstack    

    from ibvpy.mats.mats2D.mats2D_explorer_bcond import BCDofProportional
    from os.path import join
    
    from promod.simdb import SimDB
    simdb = SimDB()

    test_file = join( simdb.exdata_dir, 'tensile_tests', 'TT-9u',
                              'TT06-9u-V1.DAT' )

    ex_run = ExRun( test_file )
    
    #max_strain = 0.2*9.7599/ 1000.  #1.96666667e-03
    max_strain = ex_run.ex_type.eps_fit[-1] * 1.1
    n_steps = 4
    
    ec = {
          # overload the default configuration
          'bcond_list'  : [ BCDofProportional( max_strain = 1.0, alpha_rad = 0.0 ) ],
          'rtrace_list' : [
               RTraceGraph(name = 'stress - strain',
                           var_x = 'eps_app', idx_x = 0,
                           var_y = 'sig_app', idx_y = 0,
                           update_on = 'iteration' ),
#                           update_on = 'update' ),
                        ],
            'tline' : TLine( step = max_strain / n_steps, max = max_strain )
          }

    mats_eval = MATS2DMicroplaneDamage( 
                                   E = 24000, 
                                   nu = 0.25,                                       
                                   n_mp = 30,
                                   elastic_debug = False,
                                   stress_state = 'plane_stress',
                                   symmetrization = 'sum-type',
                                   model_version = 'compliance',
                                   phi_fn = PhiFnGeneral,
                                   )

    fitter = MATSDamageFitter \
            (
              ex_run = ex_run, 
              dim = MATS2DExplore( 
                                  mats_eval = mats_eval,
                                  explorer_config = ec,  
                  ),
                       )
    #---------------------------
    # Default settings:
    #---------------------------
    # plot the target data before running the fitter if set to 'True':
    plot_target = False
    
    if plot_target:
        import matplotlib.pyplot as plt
        plt.plot(fitter.xdata_target, fitter.ydata_target)
        plt.title('target_data')
        plt.xlabel('xdata')
        plt.ylabel('ydata')
        plt.show()

    #---------------------------
    # run fitter:
    #---------------------------
    fitter.init()
    fit_response = fitter.fit_response()

    file_name = join( simdb.matdata_dir, ex_run.ex_type.reinforcement_layout + '.mats' )
    file = open( file_name, 'w' )
    import pickle
    pickle.dump( mats_eval.phi_fn.mfn, file )
    file.close()

    #fitter.configure_traits()
    
    #---------------------------
    # basic testing of fitter methods:
    #---------------------------
    
    # set to True for basic testing of the methods:
    basic_tests = False
    
    if basic_tests:        
        fitter.run_through()
        #    fitter.tloop.rtrace_mngr.rtrace_bound_list[0].configure_traits()
        fitter.tloop.rtrace_mngr.rtrace_bound_list[0].redraw()
        last_strain_run_through = fitter.tloop.rtrace_mngr.rtrace_bound_list[0].trace.xdata[:]
        last_stress_run_through = fitter.tloop.rtrace_mngr.rtrace_bound_list[0].trace.ydata[:]
        print 'last strain (run-through) value' , last_strain_run_through
        print 'last stress (run-through) value' , last_stress_run_through
    
        fitter.tloop.reset()
        fitter.run_step_by_step()
        #fitter.tloop.rtrace_mngr.rtrace_bound_list[0].configure_traits()
        fitter.tloop.rtrace_mngr.rtrace_bound_list[0].redraw()
        last_strain_step_by_step = fitter.tloop.rtrace_mngr.rtrace_bound_list[0].trace.xdata[:]
        last_stress_step_by_step = fitter.tloop.rtrace_mngr.rtrace_bound_list[0].trace.ydata[:]
        print 'last stress (step-by-step) value', last_stress_step_by_step
        
        fitter.run_trial_step()
        fitter.run_trial_step()
        fitter.tloop.rtrace_mngr.rtrace_bound_list[0].redraw()
        strain_after_trial_steps = fitter.tloop.rtrace_mngr.rtrace_bound_list[0].trace.xdata[:]
        stress_after_trial_steps = fitter.tloop.rtrace_mngr.rtrace_bound_list[0].trace.ydata[:]
        print 'stress after trial', stress_after_trial_steps
    
        fitter.init()
        #fitter.mats2D_eval.configure_traits()
        lof = fitter.get_lack_of_fit( 1.0 )
        print '1',lof
        lof = fitter.get_lack_of_fit( 0.9 )
        print '2',lof
    
        #fitter.tloop.rtrace_mngr.configure_traits()
        fitter.run_trial_step()
    
    else:
        from ibvpy.plugins.ibvpy_app import IBVPyApp
        ibvpy_app = IBVPyApp( ibv_resource = fitter )
        ibvpy_app.main()       

if __name__ == '__main__':
    run()