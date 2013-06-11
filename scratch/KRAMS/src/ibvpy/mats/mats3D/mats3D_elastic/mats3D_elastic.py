
from enthought.traits.api import \
     Array, Bool, Callable, Enum, Float, HasTraits, \
     Instance, Int, Trait, Range, HasTraits, on_trait_change, Event, \
     implements, Dict, Property, cached_property, Delegate

from enthought.traits.ui.api import \
     Item, View, HSplit, VSplit, VGroup, Group, Spring

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

from ibvpy.mats.mats_eval import \
    IMATSEval
from ibvpy.mats.mats3D.mats3D_eval import MATS3DEval
from ibvpy.api import \
    RTrace, RTraceGraph, RTraceArraySnapshot

#---------------------------------------------------------------------------
# Material time-step-evaluator for Scalar-Damage-Model
#---------------------------------------------------------------------------

class MATS3DElastic( MATS3DEval ):
    '''
    Elastic Model.
    '''

    implements( IMATSEval )

    #---------------------------------------------------------------------------
    # Parameters of the numerical algorithm (integration)
    #---------------------------------------------------------------------------

    #---------------------------------------------------------------------------
    # Material parameters 
    #---------------------------------------------------------------------------

    E = Float( 1., #34e+3,
                 label = "E",
                 desc = "Young's Modulus",
                 auto_set = False )
    nu = Float( 0.2,
                 label = 'nu',
                 desc = "Poison's ratio",
                 auto_set = False )

    D_el = Property( Array( float ), depends_on = 'E, nu' )
    @cached_property
    def _get_D_el( self ):
        D_mtx = zeros( ( 6, 6 ), dtype = 'float_' )
        D_mtx[0, 0] = self.E / ( 1 + self.nu ) + self.E * self.nu / ( 1 + self.nu ) / ( 1 - 2 * self.nu );
        D_mtx[0, 1] = self.E * self.nu / ( 1 + self.nu ) / ( 1 - 2 * self.nu );
        D_mtx[0, 2] = self.E * self.nu / ( 1 + self.nu ) / ( 1 - 2 * self.nu );
        D_mtx[1, 0] = self.E * self.nu / ( 1 + self.nu ) / ( 1 - 2 * self.nu );
        D_mtx[1, 1] = self.E / ( 1 + self.nu ) + self.E * self.nu / ( 1 + self.nu ) / ( 1 - 2 * self.nu );
        D_mtx[1, 2] = self.E * self.nu / ( 1 + self.nu ) / ( 1 - 2 * self.nu );
        D_mtx[2, 0] = self.E * self.nu / ( 1 + self.nu ) / ( 1 - 2 * self.nu );
        D_mtx[2, 1] = self.E * self.nu / ( 1 + self.nu ) / ( 1 - 2 * self.nu );
        D_mtx[2, 2] = self.E / ( 1 + self.nu ) + self.E * self.nu / ( 1 + self.nu ) / ( 1 - 2 * self.nu );
        D_mtx[3, 3] = ( self.E / ( 1 + self.nu ) ) / 2.0
        D_mtx[4, 4] = ( self.E / ( 1 + self.nu ) ) / 2.0
        D_mtx[5, 5] = ( self.E / ( 1 + self.nu ) ) / 2.0
        return D_mtx

    # This event can be used by the clients to trigger an action upon
    # the completed reconfiguration of the material model
    #
    changed = Event

    #---------------------------------------------------------------------------------------------
    # View specification
    #---------------------------------------------------------------------------------------------

    view_traits = View( VSplit( Group( Item( 'E' ),
                                      Item( 'nu' ), ),
                                ),
                        resizable = True
                        )

    #-----------------------------------------------------------------------------------------------
    # Private initialization methods
    #-----------------------------------------------------------------------------------------------



    #-----------------------------------------------------------------------------------------------
    # Setup for computation within a supplied spatial context
    #-----------------------------------------------------------------------------------------------



    def new_cntl_var( self ):
        return zeros( 6, float_ )

    def new_resp_var( self ):
        return zeros( 6, float_ )


    #-----------------------------------------------------------------------------------------------
    # Evaluation - get the corrector and predictor
    #-----------------------------------------------------------------------------------------------

    def get_corr_pred( self, sctx, eps_app_eng, d_eps, tn, tn1 ):
        '''
        Corrector predictor computation.
        @param eps_app_eng input variable - engineering strain
        '''
        sigma = dot( self.D_el, eps_app_eng )
        return  sigma, self.D_el

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
    def _rte_dict_default( self ):
        return {'sig_app' : self.get_sig_app,
                'eps_app'  : self.get_eps_app,
                'max_principle_sig' : self.get_max_principle_sig,
                'strain_energy' : self.get_strain_energy                 }
