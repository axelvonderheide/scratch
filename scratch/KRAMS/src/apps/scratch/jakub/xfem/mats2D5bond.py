
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
from mathkit.mfn.mfn_line.mfn_line import MFnLineArray

#---------------------------------------------------------------------------
# Material time-step-evaluator for Scalar-Damage-Model
#---------------------------------------------------------------------------

class MATS2D5Bond( MATSEval ):
    '''
    Scalar Damage Model.
    '''

    implements( IMATSEval )

    E_f   = Float( 1.,#34e+3,
                 label = "E",
                 desc = "Young's Modulus Reinforcement",
                 auto_set = False )
    nu_f  = Float( 0.,#34e+3,
                 label = "nu_f",
                 desc = "Poison's ratio Reinforcement",
                 auto_set = False )
    E_m   = Float( 1.,#34e+3,
                 label = "E",
                 desc = "Young's Modulus Matrix",
                 auto_set = False )
    nu_m   = Float( 0.,#34e+3,
                 label = "nu_m",
                 desc = "Poison's ratio Matrix",
                 auto_set = False )
    G   = Float( 1.,#34e+3,
                 label = "G",
                 desc = "Shear Stiffness",
                 auto_set = False )
#    bond_fn = Trait(MFnLineArray(ydata=[0,1]),
#                    label = "Bond",
#                    desc = "Bond Function",
#                    auto_set = False)
    
    stress_state  = Enum("plane_stress","plane_strain")
    
    D_el = Property(Array(float), depends_on = 'E_f, nu_f, E_m, nu_m, stress_state')
    @cached_property
    def _get_D_el(self):
        if self.stress_state == "plane_stress":
            return self._get_D_plane_stress()
        else:
            return self._get_D_plane_strain()
    
    
    def _bond_fn_default(self):
        return MFnLineArray(xdata = [0.,1.], ydata = [0.,1.])

   # This event can be used by the clients to trigger an action upon
    # the completed reconfiguration of the material model
    #
    changed = Event                     

    #---------------------------------------------------------------------------------------------
    # View specification
    #---------------------------------------------------------------------------------------------

    #-----------------------------------------------------------------------------------------------
    # Private initialization methods
    #-----------------------------------------------------------------------------------------------
 
    #-----------------------------------------------------------------------------------------------
    # Setup for computation within a supplied spatial context
    #-----------------------------------------------------------------------------------------------

    def setup( self, sctx ):
        '''
        Intialize state variables.
        '''
        self.update_state_on = False
        
        
    def new_cntl_var(self):
        return zeros( 11, float_ )

    def new_resp_var(self):
        return zeros( 11, float_ )
        
    #-----------------------------------------------------------------------------------------------
    # Evaluation - get the corrector and predictor
    #-----------------------------------------------------------------------------------------------

    def get_corr_pred( self, sctx, eps_app_eng, d_eps, tn, tn1 ):
        '''
        Corrector predictor computation.
        @param eps_app_eng input variable - engineering strain
        '''
        #print "eps ", eps_app_eng
        
        # You print the stress you just computed and the value of the apparent E

        sigma = dot(self.D_el,eps_app_eng)
        #print "sig ", sigma
        return  sigma, self.D_el
 
    #---------------------------------------------------------------------------------------------
    # Subsidiary methods realizing configurable features
    #---------------------------------------------------------------------------------------------

    def _get_D_plane_stress( self ):
        E_s = self.E_f+self.E_m
        nu_s = (self.nu_f+self.nu_f)/2.
        E_m = self.E_m
        nu_m = self.nu_m
        E_f = self.E_f
        nu_f = self.nu_f
        G = self.G
        D_stress = zeros([11,11])
        D_stress[0][0] = E_s/(1.0-nu_s*nu_s)
        D_stress[0][1] = E_s/(1.0-nu_s*nu_s)*nu_s
        D_stress[1][0] = E_s/(1.0-nu_s*nu_s)*nu_s
        D_stress[1][1] = E_s/(1.0-nu_s*nu_s)
        D_stress[2][2] = E_s/(1.0-nu_s*nu_s)*(1.0/2.0-nu_s/2.0)
        
        D_stress[3][3] = E_m/(1.0-nu_m*nu_m)
        D_stress[3][4] = E_m/(1.0-nu_m*nu_m)*nu_m
        D_stress[4][3] = E_m/(1.0-nu_m*nu_m)*nu_m
        D_stress[4][4] = E_m/(1.0-nu_m*nu_m)
        D_stress[5][5] = E_m/(1.0-nu_m*nu_m)*(1.0/2.0-nu_m/2.0)
        
        D_stress[6][6] = E_f/(1.0-nu_f*nu_f)
        D_stress[6][7] = E_f/(1.0-nu_f*nu_f)*nu_f
        D_stress[7][6] = E_f/(1.0-nu_f*nu_f)*nu_f
        D_stress[7][7] = E_f/(1.0-nu_f*nu_f)
        D_stress[8][8] = E_f/(1.0-nu_f*nu_f)*(1.0/2.0-nu_f/2.0)
        
        D_stress[9][9]   = G
        D_stress[10][10] = G
        return D_stress

    def _get_D_plane_strain( self ):
        #TODO: adapt
        E = self.E
        nu = self.nu
        D_strain = zeros([3,3])
        D_strain[0][0] = E*(1.0-nu)/(1.0+nu)/(1.0-2.0*nu)
        D_strain[0][1] = E/(1.0+nu)/(1.0-2.0*nu)*nu
        D_strain[1][0] = E/(1.0+nu)/(1.0-2.0*nu)*nu
        D_strain[1][1] = E*(1.0-nu)/(1.0+nu)/(1.0-2.0*nu)
        D_strain[2][2] = E*(1.0-nu)/(1.0+nu)/(2.0-2.0*nu)
        return D_strain

    def get_sig_app( self, sctx, eps_app_eng ):
        # @TODO
        # the stress calculation is performed twice - it might be
        # cached.
        sig_eng, D_mtx = self.get_corr_pred( sctx, eps_app_eng, 0, 0, 0 )
        return sig_eng[:9]
    
    def get_shear( self, sctx, eps_app_eng ):
        # @TODO
        # the stress calculation is performed twice - it might be
        # cached.
        sig_eng, D_mtx = self.get_corr_pred( sctx, eps_app_eng, 0, 0, 0 )
        return sig_eng[9:]


    def get_eps_app( self, sctx, eps_app_eng ):
        return eps_app_eng[:9]
    
    def get_slip( self, sctx, eps_app_eng ):
        return eps_app_eng[9:]
    
    #---------------------------------------------------------------------------------------------
    # Response trace evaluators
    #---------------------------------------------------------------------------------------------

 
    # Declare and fill-in the rte_dict - it is used by the clients to
    # assemble all the available time-steppers.
    #
    rte_dict = Trait( Dict )
    def _rte_dict_default(self):
        return {'sig_app'      : self.get_sig_app,
                'eps_app'      : self.get_eps_app,
                'shear'        : self.get_shear,
                'slip'         : self.get_slip}
