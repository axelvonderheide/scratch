
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

class MATS1DElasticTF( MATSEval ):
    '''
    Scalar Damage Model.
    '''

    implements( IMATSEval )
     
    E_m   = Float( 1.,#34e+3,
                 label = "E_m",
                 desc = "Young's Modulus",
                 auto_set = False )

    E_f   = Float( 1.,#34e+3,
                 label = "E_f",
                 desc = "Young's Modulus",
                 auto_set = False )

    G   = Float( 1.,#34e+3,
                 label = "G",
                 desc = "Shear Modulus",
                 auto_set = False )


    # This event can be used by the clients to trigger an action upon
    # the completed reconfiguration of the material model
    #
    changed = Event                     

    #---------------------------------------------------------------------------------------------
    # View specification
    #---------------------------------------------------------------------------------------------

    view_traits = View( Item('E_m'),
                        )

 
    #-----------------------------------------------------------------------------------------------
    # Setup for computation within a supplied spatial context
    #-----------------------------------------------------------------------------------------------

    def setup( self, sctx ):
        '''
        Intialize state variables.
        '''
        self.update_state_on = False
        
        
    def new_cntl_var(self):
        return zeros( 1, float_ )

    def new_resp_var(self):
        return zeros( 1, float_ )
        
    #-----------------------------------------------------------------------------------------------
    # Evaluation - get the corrector and predictor
    #-----------------------------------------------------------------------------------------------

    def get_corr_pred( self, sctx, eps_app_eng, d_eps, tn, tn1 ):
        '''
        Corrector predictor computation.
        @param eps_app_eng input variable - engineering strain
        '''
        D_el = zeros((3,3))
        D_el[0,0] = self.E_m
        D_el[1,1] = self.E_f
        D_el[2,2] = self.G
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
#--------------------------------------------------------------------------------
# Example 
#--------------------------------------------------------------------------------

from ibvpy.core.tloop import TLoop, TLine
from ibvpy.api import BCDof
from ibvpy.core.ibvp_solve import IBVPSolve as IS

  
if __name__ == '__main__':
    # tseval for a material model
    #
    tseval  = MATS1DElastic( )
    ts = TS( tse = tseval,
             bcond_list = [ BCDof(var='u', dof = 0, value = 1.) 
                         ],
             rtrace_list = [ RTraceGraph(name = 'strain 0 - stress 0',
                                  var_x = 'eps_app', idx_x = 0,
                                  var_y = 'sig_app', idx_y = 0,
                                  update_on = 'update' )
                         ]
                         )

    # Put the time-stepper into the time-loop
    #
  
    tmax = 4.0
        # tmax = 0.0006
    n_steps = 100

    tl = TLoop( tstepper = ts,
             DT=tmax/n_steps, KMAX = 100, RESETMAX = 0,
             tline = TLine( min = 0.0,  max = tmax ) )
    
    # Put the whole stuff into the simulation-framework to map the
    # individual pieces of definition into the user interface.
    #
    from ibvpy.plugins.ibvpy_app import IBVPyApp
    app = IBVPyApp( tloop = tl )
    app.main()