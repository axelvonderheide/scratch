
from enthought.traits.api import \
     Array, Bool, Callable, Enum, Float, HasTraits, \
     Instance, Int, Trait, Range, HasTraits, on_trait_change, Event, \
     implements, Dict, Property, cached_property, Delegate

from enthought.traits.ui.api import \
     Item, View, HSplit, VSplit, VGroup, Group, Spring

# Chaco imports
from enthought.chaco.chaco_plot_editor import \
     ChacoPlotEditor, \
     ChacoPlotItem
from enthought.enable.component_editor import \
     ComponentEditor
from enthought.chaco.tools.api import \
     PanTool, SimpleZoom
from enthought.chaco.api import \
     Plot, AbstractPlotData, ArrayPlotData

#from dacwt import DAC

from numpy import \
     array, ones, zeros, outer, inner, transpose, dot, frompyfunc, \
     fabs, sqrt, linspace, vdot, identity, tensordot, \
     sin as nsin, meshgrid, float_, ix_, diag,\
     vstack, hstack, sqrt as arr_sqrt

from math import pi as Pi, cos, sin, exp, sqrt as scalar_sqrt

from scipy.linalg import eig, inv

from ibvpy.core.tstepper import \
     TStepper as TS

from ibvpy.mats.mats_eval import IMATSEval, MATSEval

from ibvpy.api import RTrace, RTraceGraph, RTraceArraySnapshot

from yield_face2D import IYieldFace2D

#---------------------------------------------------------------------------
# Material time-step-evaluator for Scalar-Damage-Model
#---------------------------------------------------------------------------

class MATS2DPlastic( MATSEval ):
    '''
    Elastic Model.
    '''

    implements( IMATSEval )

    #---------------------------------------------------------------------------
    # Parameters of the numerical algorithm (integration)
    #---------------------------------------------------------------------------
  
    stress_state  = Enum("plane strain","plane stress",)
    algorithm  = Enum("cutting plane","closest point",)
   
    #---------------------------------------------------------------------------
    # Material parameters 
    #---------------------------------------------------------------------------
    yf = Instance( IYieldFace2D,
                 label = "Yield Face",
                 desc = "Yield Face Definition")
    E   = Float( 210.e+3,
                 label = "E",
                 desc = "Young's Modulus")
    nu  = Float( 0.,
                 label = 'nu',
                 desc = "Poison's ratio" )
    K_bar = Float( 0.,
                 label = 'K',
                 desc = "isotropic softening parameter")
    H_bar = Float( 0.,
                 label = 'K',
                 desc = "kinematic softening parameter")
    tolerance = Float( 1.e-8,
                 label = 'TOL',
                 desc = "tolerance of return mapping")
    
    max_iter = Int( 20,
                    label = 'Iterations',
                    desc = "maximal number of iterations")
    
    # This event can be used by the clients to trigger an action upon
    # the completed reconfiguration of the material model
    #
    changed = Event                     

    #---------------------------------------------------------------------------------------------
    # View specification
    #---------------------------------------------------------------------------------------------

    view_traits = View( VSplit( Group(Item('yf'),
                                      Item('E'),
                                      Item('nu'),
                                      Item('tolerance'),
                                      Item('max_iter')),
                                Group( Item('stress_state', style = 'custom' ),
                                       Item('algorithm', style = 'custom' ),
                                       Spring(resizable = True),
                                       label='Configuration parameters', show_border=True,
                                       ),
                                ),
                        resizable = True
                        )

    #-----------------------------------------------------------------------------------------------
    # Private initialization methods
    #-----------------------------------------------------------------------------------------------


 
    #-----------------------------------------------------------------------------------------------
    # Setup for computation within a supplied spatial context
    #-----------------------------------------------------------------------------------------------
    def get_state_array_size( self ):
        '''
        Return number of number to be stored in state array
        @param sctx:spatial context
        '''
        return 7
 
    def setup( self, sctx ):
        '''
        Intialize state variables.
        '''
        sctx.mats_state_array = zeros( 7, float_ )
        if self.stress_state == "plane stress":
            self.D_el = self._get_D_plane_stress()
        else:
            self.D_el = self._get_D_plane_strain()
        self.H_mtx = diag([self.K_bar,self.H_bar,self.H_bar,self.H_bar])
    

    def new_cntl_var(self):
        return zeros( 3, float_ )

    def new_resp_var(self):
        return zeros( 3, float_ )

        
    #-----------------------------------------------------------------------------------------------
    # Evaluation - get the corrector and predictor
    #-----------------------------------------------------------------------------------------------

    def get_corr_pred( self, sctx, eps_app_eng, d_eps, tn, tn1 ):
        '''
        Corrector predictor computation.
        @param eps_app_eng input variable - engineering strain
        '''
        self.delta_gamma = 0.
        if sctx.update_state_on:
            print "in us"
            eps_n = eps_app_eng - d_eps
            epsilon_p, q_1, q_2 = self._get_state_variables( sctx, eps_n)
            
            sctx.mats_state_array[:3] = epsilon_p
            sctx.mats_state_array[3] = q_1
            sctx.mats_state_array[4:] =  q_2   
                              
        self.diff1s = zeros([3])
        epsilon_p, q_1, q_2 = self._get_state_variables( sctx, eps_app_eng)
        diff2ss = self.yf.get_diff2ss(eps_app_eng, self.E, self.nu, sctx)#Note: the state variables are not needed here, just gamma
        Xi_mtx = inv(inv(self.D_el) +  self.delta_gamma * diff2ss *self.f_trial)
        N_mtx_denom = sqrt(dot(dot(self.diff1s,Xi_mtx),self.diff1s))
        if N_mtx_denom == 0.:
            N_mtx = zeros(3)
        else:
            N_mtx = dot(Xi_mtx,self.diff1s)/N_mtx_denom
        D_mtx = Xi_mtx - vdot(N_mtx,N_mtx)#vdot!!
        
        # You print the stress you just computed and the value of the apparent E
        print "sigma ",self.sigma
        print "D_mtx ",D_mtx
        return  self.sigma.copy(), D_mtx
 
    #---------------------------------------------------------------------------------------------
    # Subsidiary methods realizing configurable features
    #---------------------------------------------------------------------------------------------
    def _get_state_variables( self, sctx, eps_app_eng):
        epsilon_p = sctx.mats_state_array[:3]
        q_1 = sctx.mats_state_array[3]
        q_2 = sctx.mats_state_array[4:]
       
        self.sigma =  dot(self.D_el, (eps_app_eng-epsilon_p))
        xi_trial = self.sigma - q_2
        self.f_trial = self.yf.get_f_trial(xi_trial, q_1)
        print "f_trial", self.f_trial
        int_count = 1
        while self.f_trial > self.tolerance: 
            if int_count > self.max_iter:
                print "Maximal number of iteration reached"
                break
            self.diff1s = self.yf.get_diff1s(eps_app_eng, self.E, self.nu, sctx)
            diff1q = self.yf.get_diff1q(eps_app_eng, self.E, self.nu, sctx)

            if self.stress_state == "plane stress":
                raise NotImplementedError
            else:
                if self.algorithm == "cutting plane":
                    delta_gamma_2 = self.f_trial/(dot(dot(self.diff1s,self.D_el),\
                            self.diff1s) + dot(dot(diff1q,self.H_mtx),diff1q) )
                elif self.algorithm == "closest point":
                    raise NotImplementedError
                else:
                    raise NotImplementedError
            
            self.delta_gamma += delta_gamma_2             
            epsilon_p += delta_gamma_2 * self.diff1s
                  
            q_1 += delta_gamma_2 * self.K_bar * diff1q[0]
            q_2 += delta_gamma_2 * self.H_bar * diff1q[1:]

           
            self.sigma =  dot(self.D_el,(eps_app_eng - epsilon_p))
            xi_trial = self.sigma - q_2

            self.f_trial = self.yf.get_f_trial(xi_trial, q_1)
            print "f_trial_after", self.f_trial
            int_count +=  1
        return epsilon_p, q_1, q_2

    def _get_D_plane_stress( self ):
        E = self.E
        nu = self.nu
        D_stress = zeros([3,3])
        D_stress[0][0] = E/(1.0-nu*nu)
        D_stress[0][1] = E/(1.0-nu*nu)*nu
        D_stress[1][0] = E/(1.0-nu*nu)*nu
        D_stress[1][1] = E/(1.0-nu*nu)
        D_stress[2][2] = E/(1.0-nu*nu)*(1.0/2.0-nu/2.0)
        return D_stress

    def _get_D_plane_strain( self ):
        E = self.E
        nu = self.nu
        D_strain = zeros([3,3])
        D_strain[0][0] = E*(1.0-nu)/(1.0+nu)/(1.0-2.0*nu)
        D_strain[0][1] = E/(1.0+nu)/(1.0-2.0*nu)*nu
        D_strain[1][0] = E/(1.0+nu)/(1.0-2.0*nu)*nu
        D_strain[1][1] = E*(1.0-nu)/(1.0+nu)/(1.0-2.0*nu)
        D_strain[2][2] = E*(1.0-nu)/(1.0+nu)/(2.0-2.0*nu)
        return D_strain


    #---------------------------------------------------------------------------------------------
    # Response trace evaluators
    #---------------------------------------------------------------------------------------------

    def get_sig_norm( self, sctx, eps_app_eng ):
        sig_eng, D_mtx = self.get_corr_pred( sctx, eps_app_eng, 0, 0, 0 )
        return array( [ scalar_sqrt( sig_eng[0]**2 + sig_eng[1]**2 ) ] )
    

    # Declare and fill-in the rte_dict - it is used by the clients to
    # assemble all the available time-steppers.
    #
    rte_dict = Trait( Dict )
    def _rte_dict_default(self):
        return {'sig_app'  : self.get_sig_app,
                'eps_app'  : self.get_eps_app,
                'sig_norm' : self.get_sig_norm}

if __name__ == '__main__':
    #--------------------------------------------------------------------------------
    # Example using the mats2d_explore 
    #--------------------------------------------------------------------------------
    from ibvpy.mats.mats2D.mats2D_explore import MATS2DExplore
    from yield_face2D import J2
    mats2D_explore = \
    MATS2DExplore( mats2D_eval = MATS2DPlastic(yf = J2()),
                    rtrace_list = [ RTraceGraph(name = 'strain 0 - stress 0',
                                                   var_x = 'eps_app', idx_x = 0,
                                                   var_y = 'sig_app', idx_y = 0,
                                                   record_on = 'update' ),
                                    RTraceGraph(name = 'strain 0 - strain 1',
                                                   var_x = 'eps_app', idx_x = 0,
                                                   var_y = 'eps_app', idx_y = 1,
                                                   record_on = 'update' ),
                                    RTraceGraph(name = 'stress 0 - stress 1',
                                                   var_x = 'sig_app', idx_x = 0,
                                                   var_y = 'sig_app', idx_y = 1,
                                                   record_on = 'update' ),
                                    RTraceGraph(name = 'time - sig_norm',
                                                   var_x = 'time', idx_x = 0,
                                                   var_y = 'sig_norm', idx_y = 0,
                                                   record_on = 'update' )

                                    ])
        
    mats2D_explore.tloop.eval()
    from ibvpy.plugins.ibvpy_app import IBVPyApp
    ibvpy_app = IBVPyApp( ibv_resource = mats2D_explore )
    ibvpy_app.main()
    
