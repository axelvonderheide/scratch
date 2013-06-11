
from enthought.traits.api import \
     Array, Bool, Callable, Enum, Float, HasTraits, \
     Instance, Int, Trait, Range, HasTraits, on_trait_change, Event, \
     implements, Dict, Property, cached_property, Delegate, Button, \
     Interface

from enthought.traits.ui.api import \
     Item, View, HSplit, VSplit, VGroup, Group, Spring

# Chaco imports
from enthought.chaco.chaco_plot_editor import \
     ChacoPlotEditor, \
     ChacoPlotItem
from mathkit.mfn.mfn_line.mfn_line_editor import MFnWTPlotItem
from enthought.enable.component_editor import \
     ComponentEditor
from enthought.chaco.tools.api import \
     PanTool, SimpleZoom
from enthought.chaco.api import \
     Plot, AbstractPlotData, ArrayPlotData
import enthought.units as units
from enthought.units.angle import degree, radian

#from dacwt import DAC

from numpy import \
     array, ones, zeros, outer, inner, transpose, dot, frompyfunc, \
     fabs, sqrt, linspace, vdot, identity, tensordot, \
     sin as nsin, meshgrid, float_, ix_, \
     vstack, hstack, sqrt as arr_sqrt

from math import pi as Pi, cos, sin, exp

from scipy.linalg import eig, inv

from mathkit.mfn.mfn_line.mfn_line import MFnLineArray

from ibvpy.core.tstepper import \
     TStepper as TS

from ibvpy.mats.mats_eval import IMATSEval, MATSEval

from api import RTrace, RTraceGraph

plot_item = MFnWTPlotItem("xdata", "ydata", 
                          #type_trait="plot_type",
                          type="line",

                          # Basic axis and label properties
                          show_label=False,
                          resizable=True,
                          orientation="h",
                          x_label = "strain",
                          y_label = "1 - damage",
                          y_auto = False,
                          title = 'Softening law for a microplane',
                             
                          # Plot properties
                          color = "black",
                          bgcolor = "white",
                          
                          # Border, padding properties
                          border_visible=False,
                          border_width=0,
                          padding_bg_color = "white")

#--------------------------------------------------------------------------------------
# Damage function for MDM
#--------------------------------------------------------------------------------------
class PhiFnBase( HasTraits ):
    '''
    Damage function.
    '''
    xdata = Array
    ydata = Array

    def refresh_plot(self):
        x_min, x_max = self.get_plot_range()
        x = linspace(x_min, x_max, 100 )
        phi_fn = frompyfunc( self, 2, 2 )
        y_, x_ = phi_fn( x, x )
        y = array([ v for v in y_ ], dtype='float' )
        self.xdata = x
        self.ydata = y

    def get_plot_range( self ):
        '''
        Get the range of the plot domain
        '''
        raise NotImplementedError

    def __call__( self, e_equiv, e_max, *c_list ):
        '''
        Evaluate the damage for a particular microplane.
        '''
        raise NotImplementedError

    # Default TraitsUI view
    traits_view = View( Group( Item('Epp'),
                               Item('Efp'),
                               plot_item,
                               label = 'Damage law',
                               show_border = True
                               ),
                        resizable=True,
                        width=800, height=800)

#--------------------------------------------------------------------------------------
# Damage function for MDM
#--------------------------------------------------------------------------------------
class PhiFnStrainSoftening( PhiFnBase ):
    '''
    Damage function.
    '''
    Epp = Float(desc = 'strain at the onset of damage')
    Efp = Float(desc = 'strain at total damaged')

    def get_plot_range( self ):
        return 0, self.Efp*10.

    @on_trait_change('Epp,Efp')
    def refresh_plot(self):
        super( PhiFnStrainSoftening, self ).refresh_plot()

    def __call__( self, e_equiv, e_max, *c_list ):
        '''
        Evaluate the damage for a particular microplane.
        '''
        if len( c_list ) == 0:
            c_list = [1.,1.]
        Epp = self.Epp * c_list[0]
        Efp = self.Efp * c_list[1]

        if e_equiv >= e_max:
            e_max = e_equiv
        #
        if e_max <= Epp:
            return 1.0, e_max
        else:
            return sqrt( Epp / e_max * exp( -(e_max-Epp)/Efp )), e_max


class IMPArray( Interface ):
    '''
    Microplane arrays
    '''
#--------------------------------------------------------------------------------------
# Microplane array base
#--------------------------------------------------------------------------------------
class MPArrayBase( HasTraits ):
    '''
    Manager of the microplane arrays.

    This class is responsible for the generation and initialization
    and state management of an array of microplanes. Additionally, it
    can perform the setup of damage function parameters using the
    value of the microplane integrator object.
    '''
    implements( IMPArray )
    
    n_mp = Range( 0, 6, 6,
                  label = 'Number of microplanes',
                  auto_set = False)

    # list of angles
    #
    alpha_list = Property( Array, depends_on = 'n_mp' )
    @cached_property
    def _get_alpha_list( self ):
        return array([ Pi / self.n_mp * (i - 0.5) for i in range(1,self.n_mp+1) ])


    def get_phi_list( self, e_equiv_list, sctx ):
        carr_list = self.carr_list
        n_arr = 2 + len(carr_list)
        phi_fn = frompyfunc( self.phi_fn, n_arr, 2 )
    
        #
        # damage parameter and maximum achieved equivalent strain for each
        # microplane
        # print 'phi_list, e_max_list =',phi_fn( e_equiv_list, sctx.state_array,*carr_list)
        return phi_fn( e_equiv_list, sctx.state_array,*carr_list) 


    def get_coeff_array( self ):
        raise NotImplementedError

#--------------------------------------------------------------------------------------
# Microplane Array implementation with fracture energy based damage function
#--------------------------------------------------------------------------------------
class MPArrayMDM( MPArrayBase ):
    '''
    Management of an array of microplanes.
    '''
    E   = Float( 34e+3,
                 label = "E",
                 desc = "Young's Modulus",
                 auto_set = False )
    f_t = Float( 2.8968,
                 label = 'f_t',
                 desc = 'tensile strength',
                 auto_set = False )
    G_f = Float( 0.001117 ,
                 label = 'G_f',
                 desc = 'fracture energy',
                 auto_set = False )
    md  = Float( 0.0,
                 label = 'md',
                 desc = 'factor affecting the compresive strength (explain more precisely)',
                 auto_set = False)
    h   = Float( 1.0,
                 label = 'h',
                 desc = 'element size to norm the fracture energy',
                 auto_set = False)

    #--
    # Damage function specification
    #--
    phi_fn = Instance( PhiFnStrainSoftening )
    def _phi_fn_default( self ):
        Epp, Efp = self._get_phi_fn_params()
        return PhiFnStrainSoftening( Epp = Epp, Efp = Efp )

    @on_trait_change('E,f_t,G_f,md,h')
    def update_damage_fn( self ):
        Epp, Efp = self._get_phi_fn_params()
        self.phi_fn.Epp = Epp
        self.phi_fn.Efp = Efp

    def _get_phi_fn_params( self ):
        '''
        Calculate the parameters of the damage function
        '''
        E = self.E
        f_t = self.f_t
        G_f = self.G_f
        md  = self.md
        h = self.h

        self.gamma = ( E * G_f) / ( h * f_t**2 )
        Epp = f_t / (( E * (1-md)**2 ) *( 1.95 - 0.95 / (self.gamma-1)**(0.5)))
        Efp = (G_f / ( (1-md)* h * E * Epp) + (2.13 - 1.13*md) * Epp) / (2.73-md) - Epp
        return Epp, Efp

    #---
    # Angle-dependent profiles
    #---
    fn_Epp = Instance( MFnLineArray )
    def _fn_Epp_default( self ):
        return MFnLineArray( xdata = [ 0.0, Pi  ],
                              ydata = [ 1.0, 1.0 ] )

    # coefficients
    Epp_arr = Property( Array, depends_on = 'fn_Epp,alpha_list' )
    @cached_property
    def _get_Epp_arr( self ):
        return array([ self.fn_Epp.get_value(alpha) for alpha in self.alpha_list ])

    # @todo: alex: Is the default settings for 'fn_Epf' correct?
    fn_Efp = Instance( MFnLineArray )
    def _fn_Efp_default( self ):
        return MFnLineArray( xdata = [ 0.0, Pi/2., Pi ],
                              ydata = [ 1.0, 10.0, 1.0 ] )
    
    Efp_arr = Property( Array, depends_on = 'fn_Efp,alpha_list' )
    @cached_property
    def _get_Efp_arr( self ):
        return array([ self.fn_Efp.get_value(alpha) for alpha in self.alpha_list ])

    #---
    # Coefficient array
    #---
    carr_list = Property( Array, depends_on = 'Efp_arr,Epp_arr' )
    @cached_property
    def _get_carr_list( self ):
        return [self.Epp_arr,self.Efp_arr]

    traits_view = View(Item('n_mp', style='custom',width=200),
                       Item('E'),
                       Item('f_t'),
                       Item('G_f'),
                       Item('md'),
                       Item('h'),
                       Item('phi_fn@'),
                       Item('fn_Epp@'),
                       Item('fn_Efp@'),
                       resizable = True,
                       #scrollable = True,
                       width = 0.6,
                       height = 0.6 )

#--------------------------------------------------------------------------------------
# Damage function with residual damage level for MDM
#--------------------------------------------------------------------------------------
class PhiFnStrainHardening( PhiFnBase ):
    '''
    Damage function.
    '''
    Epp = Float(5.9e-05, desc = 'strain at the onset of damage')
    Efp = Float(1.91e-04,desc = 'strain at totaly damaged state')
    Dfp = Float(1.0,desc = 'asymptotic damage level')

    def get_plot_range( self ):
        return 0.0, self.Efp*8.1

    @on_trait_change('Epp,Efp,Dfp')
    def refresh_plot(self):
        super( PhiFnStrainHardening, self ).refresh_plot()

    def __call__( self, e_equiv, e_max, *c_list ):
        '''
        Evaluate the damage for a particular microplane.
        '''
        if len( c_list ) == 0:
            c_list = [1.,1.,1.]
        Epp = self.Epp * c_list[0]
        Efp = self.Efp * c_list[1]
        Dfp = self.Dfp * c_list[2]

        if e_equiv >= e_max:
            e_max = e_equiv
        #
        if e_max <= Epp:
            return 1.0, e_max
        else:
            return (1-Dfp) * sqrt( Epp / e_max * exp( -(e_max-Epp)/Efp )) + Dfp, e_max

    # Default TraitsUI view
    traits_view = View( Group( Item('Epp'),
                               Item('Efp'),
                               Item('Dfp'),
                               plot_item,
                               label = 'Damage law',
                               show_border = True
                               ),
                        resizable=True,
                        width=800, height=800)

class MFnPolar180( HasTraits ):
    
    ralpha_range = Range( 0., 90., 90.,
                  label = 'Reinforcement range',
                  auto_set = False)
    ralpha = Range( 0., 180., 0.,
                  label = 'Reinforcement direction',
                  auto_set = False)
    salpha = Range( 0., 90., 45.,
                  label = 'Transition range',
                  auto_set = False)
    res_integ = Range(0.,1.,0.4, label = 'Residual integrity')
    fn = Instance( MFnLineArray )
    def _fn_default(self):
        xdata, ydata = self._get_fn_data()
        return MFnLineArray( xdata = xdata, ydata = ydata ) 

    @on_trait_change('ralpha_range,ralpha,salpha,res_integ')
    def _set_fn_data(self):
        self.fn.xdata, self.fn.ydata = self._get_fn_data()
    
    def _get_fn_data(self):
        ralpha       = units.convert( self.ralpha, degree, radian )
        ralpha_range = units.convert( self.ralpha_range, degree, radian )
        salpha       = units.convert( self.salpha, degree, radian )

        r_min = ralpha - ralpha_range / 2.
        r_max = ralpha + ralpha_range / 2.

        if r_min < 0:  s_min = r_min + Pi
        else:          s_min = r_min
        if r_max > Pi: s_max = r_max - Pi
        else:          s_max = r_max
        
        step_size = salpha
        rinteg = self.res_integ
        if s_min < s_max:
            self.single_interval = True
            return ([0.0, s_min-(step_size), s_min,
                    s_max, s_max+(step_size), Pi], 
                    [ 0.0, 0.0, rinteg, rinteg, 0.0, 0.0 ] )
        else:
            self.single_interval = False
            return ([0.0, s_max, s_max+(step_size),
                    s_min-(step_size), s_min, Pi],
                    [ rinteg, rinteg, 0.0, 0.0, rinteg, rinteg ])

    traits_view = View(Item('ralpha'),
                       Item('ralpha_range'),
                       Item('salpha'),
                       Item('res_integ'),
                       Item('fn@', show_label = False ),
                       resizable = True)
            
    get_value = Delegate('fn')

#--------------------------------------------------------------------------------------
# Microplane Array implementation with fracture energy based damage function
#--------------------------------------------------------------------------------------
class MPArrayCMDM( MPArrayBase ):
    '''
    Management of an array of microplanes.
    '''
    #--
    # Damage function specification
    #--
    phi_fn = Instance( PhiFnStrainHardening )
    def _phi_fn_default( self ):
        return PhiFnStrainHardening()
    #---
    # Angle-dependent profiles
    #---
    fn_Epp = Instance( MFnLineArray )
    def _fn_Epp_default( self ):
        return MFnLineArray( xdata = [ 0.0, Pi  ],
                              ydata = [ 1.0, 1.0 ] )
    Epp_arr = Property( Array, depends_on = 'fn_Epp,alpha_list' )
    @cached_property
    def _get_Epp_arr( self ):
        return array([ self.fn_Epp.get_value(alpha) for alpha in self.alpha_list ])

    fn_Efp = Instance( MFnLineArray )
    def _fn_Efp_default( self ):
        return MFnLineArray( xdata = [ 0.0, Pi/4., Pi ],
                              ydata = [ 1.0, 1.0, 1.0 ] )
    Efp_arr = Property( Array, depends_on = 'fn_Efp,alpha_list' )
    @cached_property
    def _get_Efp_arr( self ):
        return array([ self.fn_Efp.get_value(alpha) for alpha in self.alpha_list ])

    fn_Dfp = Instance( MFnPolar180 )
    def _fn_Dfp_default( self ):
        return MFnPolar180()

    Dfp_arr = Property( Array, depends_on = 'fn_Dfp.+,alpha_list' )
    @cached_property
    def _get_Dfp_arr( self ):
        return array([ self.fn_Dfp.get_value(alpha) for alpha in self.alpha_list ])

    @on_trait_change('dump_button')
    def print_data(self,event = None):
        print self.Dfp_arr

    #---
    # Coefficient array
    #---
    carr_list = Property( Array, depends_on = 'Efp_arr,Epp_arr,Dfp_arr' )
    @cached_property
    def _get_carr_list( self ):
        return [self.Epp_arr,self.Efp_arr,self.Dfp_arr]

    alpha = Range( 0., Pi, 0. )
    #--
    # Damage function specification
    #--
    phi_view = Instance( PhiFnStrainHardening )
    def _phi_view_default( self ):
        return PhiFnStrainHardening()

    @on_trait_change( 'alpha' )
    def _update_phi_view( self ):
        self.phi_view.Epp = self.phi_fn.Epp * self.fn_Epp.get_value( self.alpha )
        self.phi_view.Efp = self.phi_fn.Efp * self.fn_Efp.get_value( self.alpha )
        self.phi_view.Dfp = self.phi_fn.Dfp * self.fn_Dfp.get_value( self.alpha )

    traits_view = View(Item('n_mp', style='custom',width=200),
                       Item('alpha'),
                       #Item('phi_fn@'),
                       Item('phi_view@'),
                       #Item('fn_Epp@'),
                       #Item('fn_Efp@'),
                       Item('fn_Dfp@'),
                       #resizable = True,
                       #scrollable = True,
                       width = 0.6,
                       height = 0.6 )

    xtraits_view = View(Item('n_mp', style='custom',width=200),
                       Item('alpha'),
                       Item('phi_view@'),
                       #Item('fn_Epp@'),
                       #Item('fn_Efp@'),
                       Item('fn_Dfp@'),
                       resizable = True,
                       scrollable = True,
                       width = 0.6,
                       height = 0.6 )


if __name__ == '__main__':    

    phi_fn = PhiFnStrainHardening( Epp = 0.2, Efp = 0.6, Dfp = 0.2 )
    #phi_fn.refresh_plot()
    #phi_fn.configure_traits()

    phi_fn = MPArrayCMDM()
    phi_fn.configure_traits( view = 'xtraits_view' )

    mf = MFnLineArray( xdata = [ 0.0, Pi / 6. - 0.01, Pi / 6. + 0.01, Pi ],
                        ydata = [ 1.0, 1.0, 0.0, 0.0] )
    # mf.configure_traits( )
