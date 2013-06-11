      
from enthought.traits.api import \
     Array, Bool, Callable, Enum, Float, HasTraits, \
     Instance, Int, Trait, Range, HasTraits, on_trait_change, Event, \
     implements, Dict, Property, cached_property, Delegate, Button, \
     Interface, WeakRef, String, List, Constant, Str

# Traits UI imports
from enthought.traits.ui.api import \
     Item, View, HSplit, VSplit, VGroup, Group, Spring, TabularEditor
from enthought.traits.ui.menu \
    import OKButton, CancelButton
from enthought.traits.ui.tabular_adapter \
    import TabularAdapter

# Chaco imports
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
     
from mathkit.mfn.mfn_line.mfn_line import MFnLineArray     

from math import pi as Pi, cos, sin, exp

from scipy.linalg import eig, inv

from ibvpy.core.tstepper import \
     TStepper as TS

from ibvpy.api import RTrace, RTraceGraph

# @todo change to mfn_polar
from mathkit.mfn.mfn_polar.mfn_polar import MFnPolar

from mats2D_cmdm_phi_fn import \
    IPhiFn, PhiFnQuasiBrittle, PhiFnTensionStiffening, PhiFnGeneral

#----------------------------------------------------------------------------------
#                            Interface for Polar Damage Functions
#----------------------------------------------------------------------------------
class IPolarFn( Interface ):
    '''
    Interface to polar damage functions.
    
    They are specialized for isotropic and anisotropic specification.
    '''
#    def fit_damage_params(self):
#        '''Set the damage parameters of damage function to respect the energetic
#        characteristics of the softening process. 
#        
#        This step can involve fitting of the responce to render the prescribed 
#        fracture energy and strength for initiallly isotropic term. (parameters 
#        of the damage function are independent of the orientation). For reinforced
#        material, orientation parameters of the damage function may be involved.   
#        '''
    def get_phi_arr( self, sctx, e_max_arr ):
        '''Return a list of damage factors corresponding to the current discretization.
        '''

    def get_fracture_energy_arr( self, sctx, e_max_arr ):
        '''Integral of each individual microplane stress-strain curve.
        '''

#----------------------------------------------------------------------------------
#                                     Polar Damage Function
#----------------------------------------------------------------------------------
class PolarFnBase(HasTraits):
    '''
    Basic implementation specifies the link to the damage function.
    
    The damage function can be of different type. Therefore, a polymorphic 
    attribute is introduced here. 
    '''
    #-----------------------------------------------------------------------------------
    # Common parameters for for isotropic and anisotropic damage function specifications
    #-----------------------------------------------------------------------------------
    n_mp = Range( 0, 50, 6,
                  label = 'Number of microplanes',
                  auto_set = False)
    E   = Float( 34e+3,
                 label = "E",
                 desc = "Young's Modulus",
                 auto_set = False, enter_set = True )
    nu  = Float( 0.2,
                 label = 'nu',
                 desc = "Poison's ratio",
                 auto_set = False, enter_set = True )
    c_T = Float( 0.0,
             label = 'c_T',
             desc = 'fraction of tangential stress accounted on each microplane',
             auto_set = False, enter_set = True )

    #-----------------------------------------------------------------------------------
    # list of angles
    #-----------------------------------------------------------------------------------
    alpha_list = Property( Array, depends_on = 'n_mp' )
    @cached_property
    def _get_alpha_list( self ):
        return array([ Pi / self.n_mp * (i - 0.5) for i in range(1,self.n_mp+1) ])

    #------------------------------------------------------------------------------------
    # Damage function specification
    #------------------------------------------------------------------------------------
    phi_fn_class = Trait(  'QuasiBrittle',#General',
                          {'QuasiBrittle'      : PhiFnQuasiBrittle,
                           'General'           : PhiFnGeneral,
                           'TensionStiffening' : PhiFnTensionStiffening } )
    
    phi_fn = Property( Instance( IPhiFn ), depends_on = 'phi_fn_class' )
    @cached_property
    def _get_phi_fn(self):
        return self.phi_fn_class_()
    
    def _set_phi_fn(self, value):
        print 'settging phi_fn to',value.Epp
        if not (self.phi_fn_class_ is value.__class__):
            raise ValueError, 'class mismatch %s != %s' % ( self.phi_fn_class_, value.__class__)
        self._phi_fn = value


class IsotropicPolarFn( PolarFnBase ):
    '''
    Isotropic damage function implementation.
    
    This version does not allow anisotropic damage function specification.
    The damage evolution is however anisotropic.
    '''
    implements(IPolarFn)

    f_t = Float( 2.8968,
                 label = 'f_t',
                 desc = 'tensile strength',
                 auto_set = False, enter_set = True )
    G_f = Float( 0.001117 ,
                 label = 'G_f',
                 desc = 'fracture energy',
                 auto_set = False, enter_set = True )
    md  = Float( 0.0,
                 label = 'md',
                 desc = 'factor affecting the compresive strength (explain more precisely)',
                 auto_set = False, enter_set = True )
    h   = Float( 1.0,
                 label = 'h',
                 desc = 'element size to norm the fracture energy',
                 auto_set = False, enter_set = True )

    #-----------------------------------------------------------------------------------------------
    # Setup of damage function from the cumulative parameters
    #-----------------------------------------------------------------------------------------------
    def __init__(self, **traits ):
        super( IsotropicPolarFn, self ).__init__(**traits)
        self.phi_fn.fit_params( self.E, self.f_t, self.G_f, self.md, self.h )
        print 'in polar_fn __init__ called'
            
    @on_trait_change('phi_fn_class, E, f_t, G_f, md, h')
    def update_damage_fn( self ):
        # let the damage function to fit its damage parameters
        self.phi_fn.fit_params( self.E, self.f_t, self.G_f, self.md, self.h )
        print 'in polar_fn update_damage_fn called'

    #------------------------------------------------------------------------------------
    # Damage function specification
    #------------------------------------------------------------------------------------
    phi_fn_class = Trait( 'TensionStiffening',
                          {'QuasiBrittle'      : PhiFnQuasiBrittle,
                           'TensionStiffening' : PhiFnTensionStiffening,
                           'General'           : PhiFnGeneral        } )
    
    print 'in IsotropicPolarFn: phi_fn_class', phi_fn_class 
        
    def get_phi_arr( self, sctx, e_max_arr ):
        '''
        Return the integrity parameter for the given maximum strain
        '''
        # vectorize the damage function evaluation
        phi_fn_vect = frompyfunc( self.phi_fn.get_value, 1, 1 )
        # return the integrity parameter for each microplane
        return phi_fn_vect( e_max_arr )


    def get_fracture_energy_arr( self, sctx, e_max_arr ):
        '''
        Get the contributions of all microplanes to the macroscopic fracture energy. 
        Note: The value returned by phi_fn.get_integ must be multiplied by E to get
        the fracture energy of the considered microplane!
        '''
        integ_fn_vect = frompyfunc( self.phi_fn.get_integ, 1, 1 )
        return integ_fn_vect( e_max_arr ) * self.E
    

    traits_view = View( Item('n_mp', style='custom',width=200),
                              Item('E'),
                              Item('nu'),
                              Item('c_T'),
                              Item('f_t'),
                              Item('G_f'),
                              Item('md'),
                              Item('h'),
                              Item('phi_fn_class@'),
                              Item('phi_fn@'),
                       resizable = True,
                       #scrollable = True,
                       width = 0.6,
                       height = 0.6 )


#----------------------------------------------------------------------------------
#                                     VariedParam
#----------------------------------------------------------------------------------
class VariedParam(HasTraits):
    """
    Association between the spatial function and material parameter.
    """
    phi_fn = WeakRef( IPhiFn )

    varname = String
    
    reference_value = Property( Float )
    def _get_reference_value(self):
        return getattr( self.phi_fn, self.varname )
    
    switched_on = Bool( False )
    
    polar_fn = Instance( MFnPolar )
    @on_trait_change('switched_on')
    def reset_polar_fn(self):
        if self.switched_on:
            self.polar_fn = MFnPolar()
        else:
            self.polar_fn = None 
    
    # vectorized form of polar_fn to returning an array of coefficients
    polar_fn_vect = Property( Callable, depends_on = 'polar_fn' )
    @cached_property
    def _get_polar_fn_vect(self):
        if self.switched_on:
            return frompyfunc( self.polar_fn, 1, 1 )
        else:
            return frompyfunc( lambda angle: 1., 1, 1 )
    
    traits_view = View( Group( Item( 'varname', style = 'readonly', show_label = False ),
                               Item( 'switched_on@' ),
                               Item( 'reference_value', style = 'readonly'),
                               Item( 'polar_fn', style = 'custom', show_label=False, resizable = True ) ),
                        resizable = True,
                        height=800 )

#----------------------------------------------------------------------------------
# Tabular Adapter Definition 
#----------------------------------------------------------------------------------
class VariedParamAdapter ( TabularAdapter ):

    columns = [ ( 'Name',         'varname' ), 
                ( 'Switched on',  'switched_on' ),
                ( 'Variable',     'reference_value' ) ]
                
    font                      = 'Courier 10'
    variable_alignment        = Constant( 'right' )
    
#----------------------------------------------------------------------------------
# Tabular Editor Construction 
#----------------------------------------------------------------------------------
varpar_editor = TabularEditor(
    selected   = 'current_varpar',
    adapter    = VariedParamAdapter(),
    operations = [ 'move' ],
    auto_update = True
)

#--------------------------------------------------------------------------------------
# Microplane Array implementation with fracture energy based damage function
#--------------------------------------------------------------------------------------
class AnisotropicPolarFn( PolarFnBase ):
    '''
    Manager of the microplane arrays.

    This class is responsible for the generation and initialization
    and state management of an array of microplanes. Additionally, it
    can perform the setup of damage function parameters using the
    value of the microplane integrator object.
    '''
    implements(IPolarFn)

    # possible realizations of damage functions
    #
    phi_fn_class = Trait('TensionStiffening',
                         {'TensionStiffening' : PhiFnTensionStiffening } )

    varied_params = List( Str, [] )
    
    #-----------------------------------------------------------------------------------------------
    # Management of spatially varying parameters depending on the value of mats_eval
    #-----------------------------------------------------------------------------------------------
    varpars = Dict
    def _varpars_default(self):
        return self._get_varpars()
    
    @on_trait_change( 'phi_fn,varied_params' )
    def _update_varpars(self):
        self.varpars = self._get_varpars()
    
    def _get_varpars( self ):
        '''
        reset the varpar list according to the current phi_fn object.
        '''
        params = self.phi_fn.identify_parameters()
        varset = {}
        for key in params:
            par_val = getattr( self.phi_fn , key )
            varset[key] = VariedParam(phi_fn = self.phi_fn,
                                      varname = key)
#                                      polar_fn = MFnPolar180() )
            if key in self.varied_params:
                varset[key].switched_on = True
        return varset

    varpar_list = Property( List( VariedParam ), depends_on = 'varpars' )
    @cached_property
    def _get_varpar_list(self):
        return [self.varpars[key] 
                for key in self.phi_fn.identify_parameters() ]

    # variable selectable in the table of varied params (just for viewing) 
    current_varpar = Instance( VariedParam )
    def _current_varpar_default(self):
        return self.varpar_list[0]
    
    @on_trait_change('phi_fn')
    def set_current_varpar(self):
        self.current_varpar = self.varpar_list[0]

    #-----------------------------------------------------------------------------------------------
    # Get the damage state for all microplanes
    #-----------------------------------------------------------------------------------------------
    def get_phi_arr( self, sctx, e_max_arr ):
        '''
        Return the damage coefficients
        '''
        # gather the coefficients for parameters depending on the orientation
        carr_list = [self.varpars[key].polar_fn_vect( self.alpha_list ) 
                     for key in self.phi_fn.identify_parameters() ]  
        # vectorize the damage function evaluation
        n_arr = 1 + len(carr_list)
        phi_fn_vect = frompyfunc( self.phi_fn.get_value, n_arr, 1 )       
        # damage parameter for each microplane
        return phi_fn_vect( e_max_arr, *carr_list)
                 
                 
    def get_fracture_energy_arr( self, sctx, e_max_arr ):
        '''
        Return the fracture energy contributions
        '''
        carr_list = [self.varpars[key].polar_fn_vect( self.alpha_list ) 
                     for key in self.phi_fn.identify_parameters() ]
        # vectorize the damage function evaluation
        n_arr = 1 + len(carr_list)
        phi_fn_vect = frompyfunc( self.phi_fn.get_integ, n_arr, 1 )
        return self.E * phi_fn_vect( e_max_arr, *carr_list )
                         
                         
    traits_view = View(Group( Item('n_mp', style='custom',width=200),
                              Item('E'),
                              Item('nu'),
                              Item('c_T'),                              
                              Item('phi_fn_class@'),
                              Item('phi_fn@', show_label = False ),
                              label = 'Polar damage function'),
                        Group(
                              HSplit( 
                                     Item('varpar_list', show_label = False, editor = varpar_editor ),
                                     Item('current_varpar', show_label = False, style = 'custom', resizable = True),
                                     ),
                              label = 'Angle-dependent variations'
                              ),                              
                       resizable = True,
                       #scrollable = True,
                       width = 0.6,
                       height = 0.6 )


if __name__ == '__main__':    

    phi_fn_brittle = PhiFnQuasiBrittle()
#    phi_fn_brittle = PhiFnQuasiBrittle( Epp = 0.255, Efp = 0.6 )
#    phi_fn_ductile = PhiFnTensionStiffening( Epp = 0.2, Efp = 0.6, Dfp = 0.2 )

    phi_fn_brittle_array = IsotropicPolarFn( phi_fn_class = 'QuasiBrittle', 
                                             phi_fn = phi_fn_brittle )
#    phi_fn_brittle_array = IsotropicPolarFn()

    
    print '2'
    phi_fn_brittle_array.phi_fn = phi_fn_brittle
    print 'phi_fn.Epp',phi_fn_brittle_array.phi_fn.Epp
    phi_fn_brittle_array.phi_fn.Efp = 0.888888
    print '3'
    phi_fn_brittle_array.configure_traits()

    from math import pi
    
    phi_fn_ductile_array = AnisotropicPolarFn( varied_params = ['Dfp'] )
    phi_fn_ductile_array.varied_params = ['Dfp']
    phi_fn_ductile_array.varpars['Dfp'].polar_fn.set( phi_residual = 1.0,
                                                      phi_quasibrittle = 0.0,
                                                      delta_trans = pi/4. ,
                                                      delta_alpha = pi/4. ) 
    print phi_fn_ductile_array.varpars

    phi_fn_ductile_array.configure_traits()
