
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
     sin as nsin, meshgrid, float_, ix_, \
     vstack, hstack, sqrt as arr_sqrt

from math import pi as Pi, cos, sin, exp, sqrt as scalar_sqrt

from scipy.linalg import eig, inv

from ibvpy.core.tstepper import \
     TStepper as TS

from ibvpy.mats.mats_eval import IMATSEval, MATSEval

from ibvpy.api import RTrace, RTraceGraph, RTraceArraySnapshot


#---------------------------------------------------------------------------
# Material time-step-evaluator for Scalar-Damage-Model
#---------------------------------------------------------------------------

class MA1DPlastic( MATSEval ):
    '''
    Scalar Damage Model.
    '''

    implements( IMATSEval )

    E           = Float( 1.,#34e+3,
                         label = "E",
                         desc = "Young's Modulus",
                         auto_set = False )
      
    epsilon_p_n = Float(0.,
                        label = "eps_p",
                        desc = "plastic strain",
                        auto_set = False )
    
    K_bar       = Float(0.,
                        label = "K_bar",
                        desc = "isotropic hardening component",
                        auto_set = False )
    
    H_bar       = Float(0.,
                        label = "H_bar",
                        desc = "kinematic hardening component",
                        auto_set = False )
    
#    h           = Float([[K_bar,0.,0.,0.,0.,0.,0.],
#                         [0.,H_bar,0.,0.,0.,0.,0.],
#                         [0.,0.,H_bar,0.,0.,0.,0.],
#                         [0.,0.,0.,H_bar,0.,0.,0.],
#                         [0.,0.,0.,0.,H_bar,0.,0.],
#                         [0.,0.,0.,0.,0.,H_bar,0.],
#                         [0.,0.,0.,0.,0.,0.,H_bar]], 
#                         label = "K",
#                         desc = "plastic modulus",
#                         auto_set = False )
    
    stiffness  = Enum("secant","algorithmic")
    
    #return_mapping = Enum("Cutting-Plane Algorithm","Closest Point Projection")
    
                    
    # This event can be used by the clients to trigger an action upon
    # the completed reconfiguration of the material model
    #
    changed = Event                     

    #---------------------------------------------------------------------------------------------
    # View specification
    #---------------------------------------------------------------------------------------------

    view_traits = View( Item('E'),
                        )

    #-----------------------------------------------------------------------------------------------
    # Private initialization methods
    #-----------------------------------------------------------------------------------------------
    def __init__( self, **kwtraits ):
        '''
        Subsidiary arrays required for the integration.
        they are set up only once prior to the computation
        '''
        super( MA1DPlastic, self ).__init__( **kwtraits )
 
    #-----------------------------------------------------------------------------------------------
    # Setup for computation within a supplied spatial context
    #-----------------------------------------------------------------------------------------------

    def get_state_array_size( self ):
        '''
        Give back the nuber of floats to be saved
        @param sctx:spatial context
        '''
        return 2

    def setup( self, sctx ):
        '''
        Intialize state variables.
        '''
        self.update_state_on = False
        sctx.mats_state_array = zeros( 2, float_ )
        
    def new_cntl_var(self):
        return zeros( 1, float_ )

    def new_resp_var(self):
        return zeros( 1, float_ )
        
    #-----------------------------------------------------------------------------------------------
    # Evaluation - get the corrector and predictor
    #-----------------------------------------------------------------------------------------------

    def get_corr_pred( self, sctx, eps_app_eng, tn, tn1 ):
        '''
        Corrector predictor computation.
        @param eps_app_eng input variable - engineering strain
        '''

        E   = self.E
        D_el = array([[E]])
        sigma = dot( D_el, eps_app_eng )
        
        # You print the stress you just computed and the value of the apparent E

        return  sigma, D_el
 
    #---------------------------------------------------------------------------------------------
    # Subsidiary methods realizing configurable features
    #---------------------------------------------------------------------------------------------

    
    #---------------------------------------------------------------------------------------------
    # Response trace evaluators
    #---------------------------------------------------------------------------------------------

    # Declare and fill-in the rte_dict - it is used by the clients to
    # assemble all the available time-steppers.
    #
    rte_dict = Trait( Dict )
    def _rte_dict_default(self):
        return {'sig_app'      : self.get_sig_app,
                'eps_app'      : self.get_eps_app}

    