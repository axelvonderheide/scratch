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

from enthought.traits.api import Float
from ibvpy.mats.mats2D.mats2D_explore import MATS2DExplore
from numpy import copy

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

    step_size = Float( 0.002 / 30 )
    
    tmax = Float( 0.002 )

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
            current_time = tmax
            tmax = current_time + self.step_size
            self.tloop.tline.max = tmax
            self.tloop.eval()

    def run_n_trial_steps( self, n ):
        self.tloop.tline.step = self.step_size
        current_time = self.tloop.t_n
        print 'current time', current_time
        tmax = current_time + self.step_size
        self.tloop.tline.max = tmax
        self.U_n = copy( self.tloop.U_n )
        print 'tmax', tmax
        self.tloop.debug = False
        for i in range( n ):
            self.tloop.t_n = current_time
            self.tloop.U_n = copy( self.U_n )
            self.tloop.eval()
            self.tloop.tstepper.sctx.update_state_on = False


    def run_trial_step( self ):
        '''Run the computation one step starting from the
        current time t_n to iterate the value for phi.
        NOTE: The trial step does not update 'U_n' or 't_n'!
        '''
        self.tloop.tline.step = self.step_size
        current_time = self.tloop.t_n
        print 'trail step: current time = ', current_time
        tmax = current_time + self.step_size
        self.tloop.tline.max = tmax
        self.U_n = copy( self.tloop.U_n )
        print 'trail step: tmax = ', tmax
        self.tloop.debug = False
        self.tloop.t_n = current_time
        self.tloop.U_n = copy( self.U_n )
        self.tloop.eval()
        self.tloop.tstepper.sctx.update_state_on = False

    def run_one_step(self):
        '''Run the computation one step starting from the
        current time t_n with the iterated value for phi
        NOTE: The calculated step does update 'U_n' or 't_n'!
        '''
        self.tloop.tline.step = self.step_size
        current_time = self.tloop.t_n
        print 'one step: current time = ', current_time
        tmax = current_time + self.step_size
        print 'one step: tmax = ', tmax
        self.tloop.eval()
        self.tloop.tstepper.sctx.update_state_on = True







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
    
#    def get_value( self, e_max, *c_list ):
#        '''
#        This method is used for fitting of the damage function. 
#        It uses the method 'get_value' from 'Instance(MfnLineArray)'
#        which is constructed in the base class 'PhiFnBase'. 
#        '''
#        return self.mfn.get_value( e_max )    
#    
#    e_max_value = max( e_max_arr_new )
#    
#    '''set default values: the computation starts with time=0.
#    For this step 'phi_fn.get_value' returns 1.0'''
#    x = array([0., 1.])
#    y = array([1., 1.])
#    self.polar_fn.phi_fn.mfn.set( xdata = x, ydata = y )
#    self.polar_fn.phi_fn.mfn.data_changed = True
##        print 'setup xdata for 1. step: ', x
##        print 'setup ydata for 1. step: ', y
#
#    '''If the value for 'e_max_value' changes a new pair is added 
#    in the damage function defined by 'xdata' and 'ydata'. The new point 
#    'e_max_value' consists of the new 'e_max_value' and the 
#    corresponding value for 'phi' obtained from fitting. 
#    '''
#    x_ = self.polar_fn.phi_fn.mfn.xdata
#    y_ = self.polar_fn.phi_fn.mfn.ydata
#    if len(x_) == 2 and ( x_ == array([0., 1.])).all():
#        x = hstack([ x_[0], self.e_max_value ])
#        y = hstack([ y_[0], sin(self.e_max_value) ])
#        self.polar_fn.phi_fn.mfn.set( xdata = x, ydata = y )
#        self.polar_fn.phi_fn.mfn.data_changed = True
##            print 'update xdata for 2. step: ', x
##            print 'update ydata for 2. step: ', y
#    else:
#        x = hstack([x_, self.e_max_value])
#        y = hstack([y_, sin(self.e_max_value)])
#        self.polar_fn.phi_fn.mfn.set( xdata = x, ydata = y )
#        self.polar_fn.phi_fn.mfn.data_changed = True
##            print 'update xdata after 2. step: ', x
##            print 'update ydata after 2. step: ', y
    
    
        
    mdm = MA2DMicroplaneDamage( elastic_debug = False )
    
    fitter = MATS2DDamageFitter \
            (          mats2D_eval = mdm,
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


        
#    fitter.run_through()
##    fitter.tloop.rtrace_mngr.rtrace_bound_list[1].configure_traits()
#    fitter.tloop.rtrace_mngr.rtrace_bound_list[1].redraw()
#    last_strain_run_through = fitter.tloop.rtrace_mngr.rtrace_bound_list[1].trace.xdata[:]
#    last_stress_run_through = fitter.tloop.rtrace_mngr.rtrace_bound_list[1].trace.ydata[:]

#    fitter.tloop.reset = True
#    fitter.run_step_by_step()
#    #fitter.tloop.rtrace_mngr.rtrace_bound_list[1].configure_traits()
#    fitter.tloop.rtrace_mngr.rtrace_bound_list[1].redraw()
#    last_strain_step_by_step = fitter.tloop.rtrace_mngr.rtrace_bound_list[1].trace.xdata[:]
#    last_stress_step_by_step = fitter.tloop.rtrace_mngr.rtrace_bound_list[1].trace.ydata[:]
    
    
#    fitter.run_n_trial_steps( 2 )
#    fitter.tloop.rtrace_mngr.rtrace_bound_list[1].redraw()
#    strains_after_trial_steps = fitter.tloop.rtrace_mngr.rtrace_bound_list[1].trace.xdata[:]
#    stresses_after_trial_steps = fitter.tloop.rtrace_mngr.rtrace_bound_list[1].trace.ydata[:]


    
    
    fitter.tmax      = 0.002 / 1000.
    fitter.step_size = 0.002 / 1000.
    fitter.tloop.reset = True
    fitter.run_through()
    fitter.tmax      = 0.002
    fitter.step_size = 0.002 / 10.

    # fitting data (from the micro-model or from the experiment)
#    xdata_fit = array([0., 1.0])
#    ydata_fit = array([0., 1.0])

    xdata_fit = array([0.00000000e+00,   3.33333333e-05,   6.66666667e-05,   1.00000000e-04,
                       1.33333333e-04,   1.66666667e-04,   2.00000000e-04,   2.33333333e-04,
                       2.66666667e-04,   3.00000000e-04,   3.33333333e-04,   3.66666667e-04,
                       4.00000000e-04,   4.33333333e-04,   4.66666667e-04,   5.00000000e-04,
                       5.33333333e-04,   5.66666667e-04,   6.00000000e-04,   6.33333333e-04,
                       6.66666667e-04,   7.00000000e-04,   7.33333333e-04,   7.66666667e-04,
                       8.00000000e-04,   8.33333333e-04,   8.66666667e-04,   9.00000000e-04,
                       9.33333333e-04,   9.66666667e-04,   1.00000000e-03,   1.03333333e-03,
                       1.06666667e-03,   1.10000000e-03,   1.13333333e-03,   1.16666667e-03,
                       1.20000000e-03,   1.23333333e-03,   1.26666667e-03,   1.30000000e-03,
                       1.33333333e-03,   1.36666667e-03,   1.40000000e-03,   1.43333333e-03,
                       1.46666667e-03,   1.50000000e-03,   1.53333333e-03,   1.56666667e-03,
                       1.60000000e-03,   1.63333333e-03,   1.66666667e-03,   1.70000000e-03,
                       1.73333333e-03,   1.76666667e-03,   1.80000000e-03,   1.83333333e-03,
                       1.86666667e-03,   1.90000000e-03,   1.93333333e-03,   1.96666667e-03,
                       1.96666667e-03])
    
    ydata_fit = array([0.,          1.18055556,  2.32136204,  2.87776364,  3.38092992,  3.663105,
                       3.75724285,  3.84216067,  3.92147416,  3.997584  ,  4.07213255,  4.1462654,
                       4.22079085,  4.29628189,  4.37314373,  4.4516601 ,  4.53202557,  4.61436863,
                       4.69876847,  4.78526735,  4.87387975,  4.96459931,  5.0574041 ,  5.15226056,
                       5.24912666,  5.3479542 ,  5.44869074,  5.55128092,  5.65566766,  5.76179299,
                       5.86959871,  5.97902693,  6.09002044,  0.74043954,  0.74661105,  0.75304179,
                       0.75972494,  0.7666537 ,  0.77382139,  0.7812214 ,  0.78884722,  0.79669249,
                       0.80475096,  0.81301651,  0.82148317,  0.83014509,  0.83899658,  0.84803208,
                       0.85724617,  0.86663356,  0.87618913,  0.88590785,  0.89578486,  0.90581541,
                       0.91599489,  0.9263188 ,  0.93678277,  0.94738257,  0.95811405,  0.9689732,
                       0.9689732 ])



   
    # interpolation function for fitting data
    mfn_line_array_fit = MFnLineArray()
    mfn_line_array_fit.set( xdata = xdata_fit, ydata = ydata_fit )
    mfn_line_array_fit.data_changed = True
    
    
    
    # get the current strain:
    eps_app_tn1 = copy( fitter.tloop.U_k )
    print 'eps_app_tn1', eps_app_tn1    
    # get the current sctx:
    sctx = fitter.tloop.tstepper.sctx
    print 'sctx', sctx    

    # for fitting use default method 'get_value' of 'MFnLineArray' as 'phi_fn':
    fitter.mats2D_eval.polar_fn.phi_fn.mfn.get_value = MFnLineArray().get_value

    e_max_value_new = max( fitter.mats2D_eval._get_state_variables( sctx, eps_app_tn1 ) )
    print 'e_max_value_new', e_max_value_new    

    phi_value_old = float(1.0)
    phi_value_new = float(1.0)

    # setup phidata (only for eps_tn = array([0., 0., 0.])):
    # set default values: the computation starts with time = 0.
    # For this step 'phi_fn.get_value' returns 1.0
    x = array([ 0., e_max_value_new ])
    y = array([ 1.,  phi_value_old  ])
    fitter.mats2D_eval.polar_fn.phi_fn.mfn.set( xdata = x, ydata = y )
    fitter.mats2D_eval.polar_fn.phi_fn.mfn.data_changed = True
    print 'setting up xdata: ', x
    print 'setting up ydata: ', y


    # outer loop:
    for n_steps in range(10):

        # inner loop:
        KMAX = 100
        k = 1
        while k < KMAX:   
            print 'inner loop: IT:', k
            fitter.run_n_trial_steps( 1 )
            fitter.tloop.rtrace_mngr.rtrace_bound_list[1].redraw()
            # get calculated sig_app_xx based on current phi_value:
            sig_app_trial = fitter.tloop.rtrace_mngr.rtrace_bound_list[1].trace.ydata[-1]
            # get sig_app_xx to be fitted:
            sig_app_fit = mfn_line_array_fit.get_value(eps_app_tn1[0])
            
            fitting_ratio = sig_app_trial / sig_app_fit
            if fitting_ratio > 1:
                phi_value_new = phi_value_old * 0.8  
            elif fitting_ratio > 1:
                phi_value_new = phi_value_old * 1.2  
                
            # adjusting phi_value:
            x = hstack([ x[:-1], e_max_value_new ])
            y = hstack([ y[:-1],  phi_value_new  ])
            fitter.mats2D_eval.polar_fn.phi_fn.mfn.set( xdata = x, ydata = y )
            fitter.mats2D_eval.polar_fn.phi_fn.mfn.data_changed = True
            print 'adjusting phi_value to: ', phi_value_new
            if fitting_ratio < 1.01:
                fitter.run_one_step()
                break
            else:
                phi_value_old = phi_value_new
                k += 1
        
        # get the maximum value for the microplane strains based on 
        # the strains from the current step and the state variables
        e_max_value_new = max( fitter.mats2D_eval._get_state_variables( sctx, eps_app_tn1 ) )
        print 'e_max_value_new', e_max_value_new        

        # updating phidata:
        x = hstack([ x, e_max_value_new ])
        y = hstack([ y,  phi_value_new  ])
        fitter.mats2D_eval.polar_fn.phi_fn.mfn.set( xdata = x, ydata = y )
        fitter.mats2D_eval.polar_fn.phi_fn.mfn.data_changed = True
        print 'update xdata: ', x
        print 'update ydata: ', y

    

 
#    e_max_value = max( fitter.tloop.tstepper.sctx.mats_state_array )





     
     
     
#    print 'last stress (run-through) value' , last_stress_run_through
#    print 'last stress (step-by-step) value', last_stress_step_by_step
#    print 'stresses after trial'            , stresses_after_trial_steps

#    fitter.mats2D_eval.polar_fn.n_mp = 6
#
#
#    

#    fitter.run_trial_step()
    from ibvpy.plugins.ibvpy_app import IBVPyApp
    ibvpy_app = IBVPyApp( ibv_resource = fitter )
    ibvpy_app.main()       
#