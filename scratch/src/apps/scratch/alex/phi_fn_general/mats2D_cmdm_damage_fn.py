
from enthought.traits.api import \
     Array, Bool, Callable, Enum, Float, HasTraits, \
     Instance, Int, Trait, Range, HasTraits, on_trait_change, Event, \
     implements, Dict, Property, cached_property, Delegate, Button, \
     Interface, WeakRef, String, List, Constant

# Traits UI imports
from enthought.traits.ui.api import \
     Item, View, HSplit, VSplit, VGroup, Group, Spring, TabularEditor
from enthought.traits.ui.menu \
    import OKButton, CancelButton
from enthought.traits.ui.tabular_adapter \
    import TabularAdapter

from mathkit.mfn.mfn_line.mfn_matplotlib_editor import MFnMatplotlibEditor
from mathkit.mfn.mfn_line.mfn_plot_adapter import MFnPlotAdapter

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
     sin as nsin, meshgrid, float_, ix_, where, trapz,\
     vstack, hstack, sqrt as arr_sqrt
     
from mathkit.mfn.mfn_line.mfn_line import MFnLineArray     

from math import pi as Pi, cos, sin, exp

from scipy.linalg import eig, inv

from ibvpy.api import RTrace, RTraceGraph
from ibvpy.core.tstepper import TStepper as TS

#mfn_editor = MFnMatplotlibEditor( \
#            adapter = MFnPlotAdapter( label_x = 'strain',
#                                      label_y = 'integrity',
#                                      title = 'Softening law for a microplane',
#                          # Plot properties
#                          color = "black",
#                          bgcolor = "white",
#                          
#                          # Border, padding properties
#                          border_visible=False,
#                          border_width=0,
#                          padding_bg_color = "white"                                      
#                                      ) )

class IPhiFn( Interface ):
    '''Interface to adamage function representation
    '''
    def identify_parameters(self):
        '''Return the set of parameters defining the respective damage function.
        '''
    def __call__( self, e_max, *c_list ):
        '''return the value of the damage function for the,
        maximum strain achieved so far and list of coefficients for material poarameters
        '''

#--------------------------------------------------------------------------------------
# Damage function for MDM
#--------------------------------------------------------------------------------------
class PhiFnBase( HasTraits ):
    '''
    Damage function.
    '''
        
    mfn = Instance( MFnLineArray )
    def _mfn_default( self ):
        return MFnLineArray( xdata = [0,1], ydata = [1,1] )

    def __init__(self, **args ):
        super( PhiFnBase, self ).__init__(**args)
        self.refresh_plot()

    def refresh_plot(self):
        x_min, x_max = self.get_plot_range()
        x = linspace(x_min, x_max, 100 )
        phi_fn = frompyfunc( self, 1, 1 )
        x_ = x
        y_ = phi_fn( x )
        y = array([ v for v in y_ ], dtype='float' )
        self.mfn.set( xdata = x, ydata = y )
        self.mfn.data_changed = True

    def get_plot_range( self ):
        '''
        Get the range of the plot domain
        '''
        raise NotImplementedError

    def identify_parameters(self):
        '''
        Extract the traits that are of type Float 
        '''
        params = []
        for name, trait in self.traits().items():
            if trait.is_trait_type(Float):
                params.append( name )
        print 'params', params
        return params

    def fit_params( self, *params ):
        '''Possiblity to adapt the microplane-related 
        material paramters based on the integral characteric specification.
        '''
        return
    
    # Default TraitsUI view
    traits_view = View( Group( #Item('Epp'),
                               #Item('Efp'),
                               Item( 'mfn', show_label = False ),#, editor = mfn_editor ),
                               label = 'Damage law',
                               show_border = True
                               ),
                        resizable = True,
                        width=800, height=800 )

#--------------------------------------------------------------------------------------
# Piecewise Linear damage function for MDM
#--------------------------------------------------------------------------------------
class PhiFnGeneral( PhiFnBase ):

    implements( IPhiFn )

    def get_value( self, e_max, *c_list ):
        '''
        Evaluate the integrity of a particular microplane.
        '''
        return self.mfn.get_value( e_max )


    def get_integ( self, e_max, *c_list ):
        '''
        Evaluate the integrity of a particular microplane.
        '''
        # get the data that defines PhiFnGeneral
        # (fitted from experiment)
        _xdata = self.mfn.xdata
        _ydata = self.mfn.ydata

        # get the values smaller then the current e_max
        _xdata_ = _xdata[ where( _xdata < e_max )[0] ]
        _ydata_ = _ydata[ :len(_xdata_) ]
#        print '_xdata_' , _xdata_ 

        # add the value pair for e_max
        _xdata_emax = hstack([_xdata_, e_max])
        _ydata_emax = hstack([_ydata_, self.mfn.get_value(e_max)])
#        print '_xdata_emax' , _xdata_emax 
        
        # assume an uncoupled relation (e.g. direct link) between the microplane
        # strains (e) and microplane stresses (s), e.g. s = phi * E * phi * e;
        # The methode 'get_integ' returns the value without the young's modulus,
        # which is multiplied in 'PhiFnPolar';
        # _ydata_integ = phi * phi * e
        #
        # @todo: this is only approximately true!; the correct evaluation 
        # takes the version consistend (stiffness/compliance) pairs for
        # the microplane strains and stresses (work conjugates) 
        _ydata_integ = _ydata_emax * _ydata_emax * _xdata_emax
        
        # integral under the stress-strain curve
        E_t = trapz(_ydata_integ,_xdata_emax)
        # area of the stored elastic energy  
        U_t = 0.0
        if len(_xdata_emax) != 0:
            U_t = 0.5 * _ydata_integ[-1] * _xdata_emax[-1]     
#        print 'E_t', E_t
#        print 'U_t', U_t        
#        print 'E_t - U_t', E_t - U_t
        return E_t - U_t                
#        return self.mfn.integ_value

    def get_plot_range( self ):
        return self.mfn.xdata[0], self.mfn.xdata[-1]
    
    def __call__( self, e_max, *c_list ):
        return self.get_value( e_max, *c_list )

#--------------------------------------------------------------------------------------
# Damage function for MDM
#--------------------------------------------------------------------------------------
class PhiFnQuasiBrittle( PhiFnBase ):
    '''
    Damage function.
    '''

    implements( IPhiFn )
    
    # @todo: verify if the definition of a default value is necessary:
    #        otherwise a ZeroDevisionError occurred because a refresh plot
    #        is called before the parameters can be fitted   
    Epp = Float(1.0, desc = 'strain at the onset of damage', enter_set = True, auto_set = False )
    Efp = Float(1.0, desc = 'strain at total damaged', enter_set = True, auto_set = False )

    def identify_parameters(self):
        return ['Epp','Efp']

    def get_plot_range( self ):
        #return 0, 3.
        return 0, self.Epp*20.

    @on_trait_change('Epp, Efp')
    def refresh_plot(self):
        print 'refresh plot in QuasiBrittle'
        super( PhiFnQuasiBrittle, self ).refresh_plot()

    def fit_params( self, *params ):
        '''
        Calculate the parameters of the damage function
        '''
        E, f_t, G_f, md, h = params

        self.gamma = ( E * G_f) / ( h * f_t**2 )
        self.Epp = f_t / (( E * (1-md)**2 ) *( 1.95 - 0.95 / (self.gamma-1)**(0.5)))
        self.Efp = (G_f / ( (1-md)* h * E * self.Epp) + 
                    (2.13 - 1.13*md) * self.Epp) / (2.73-md) - self.Epp
        print 'fit params: self.Epp', self.Epp 
        print 'fit params: self.Efp', self.Efp 

    def get_integ( self, e_max, *c_list ):
        '''
        OBSOLETE method - was used for decoupled evaluation of fracture 
        energy contribution of the microplane.

        The method returns the value of the following integral: 
        int( Phi(e_max~)**2 * e_max~, e_max~ = 0..e_max )
        The value corresponds to the fracture energy of the considered microplane 
        divided by E. (Note: For a damage function Phi(e_max) the microplane stress-strain curve
        evaluates to s = Phi(e_max)*E*Phi(e_max)*e_max.) 
        '''
        if len( c_list ) == 0:
            c_list = [1.,1.]
        Epp = self.Epp * c_list[0]
        Efp = self.Efp * c_list[1]
        #
        # cf. derivation in Maple 'relation_between_Gf-ft_and_Epp-Efp
        # devide the result by E, i.e. the returned value does NOT include E
        if e_max <= Epp:
            return 0     
        else:
            return -0.5 * Epp * (- Epp - 2.0*Efp + 2.0*Efp * exp( ((-e_max + Epp) / Efp))) \
                   -0.5 *  e_max * Epp * exp(- ((e_max - Epp) / Efp));

        
    def get_value( self, e_max, *c_list ):
        '''Evaluate the integrity of a particular microplane 
        based on an exponential softening law.
        '''
        if len( c_list ) == 0:
            c_list = [1.,1.]
        Epp = self.Epp * c_list[0]
        Efp = self.Efp * c_list[1]
        if e_max <= Epp:
            return 1.0
        else:
            return sqrt( Epp / e_max * exp( -(e_max-Epp)/Efp ))    


    def __call__( self, e_max, *c_list ):
        return self.get_value( e_max, *c_list )

#--------------------------------------------------------------------------------------
# Damage function with residual damage level for MDM
#--------------------------------------------------------------------------------------
class PhiFnTensionStiffening( PhiFnBase ):
    '''
    Damage function.
    '''

    implements( IPhiFn )
    
    Epp = Float(5.9e-05, desc = 'microplane strain at the onset of damage',
                enter_set = True, auto_set = False )
    Efp = Float(1.91e-04,desc = 'microplane strain at totaly damaged state',
                enter_set = True, auto_set = False )
    Dfp = Float(0.3,desc = 'asymptotic damage level',
                enter_set = True, auto_set = False )
    Elimit = Float( 8.00e-02,desc = 'microplane strain at ultimate failure',
                enter_set = True, auto_set = False )

    def identify_parameters(self):
        return ['Epp','Efp','Dfp','Elimit' ]

    def get_plot_range( self ):
#        return 0.0, self.Efp*8.1
        return 0.0, self.Elimit*1.2

    @on_trait_change('Epp,Efp,Dfp,Elimit')
    def refresh_plot(self):
        super( PhiFnTensionStiffening, self ).refresh_plot()
  
    def get_integ( self, e_max, *c_list ):
        '''
        OBSOLETE method - was used for decoupled evaluation of fracture 
        energy contribution of the microplane.

        The method returns the value of the following integral: 
        int( Phi(e_max~)**2 * e_max~, e_max~ = 0..e_max )
        The value corresponds to the fracture energy of the considered microplane 
        divided by E. (Note: For a damage function Phi(e_max) the microplane stress-strain curve
        evaluates to s = Phi(e_max)*E*Phi(e_max)*e_max.) 
        '''
        if len( c_list ) == 0:
            c_list = [1.,1.]
        Epp = self.Epp * c_list[0]
        Efp = self.Efp * c_list[1]
        Dfp = self.Dfp * c_list[2]
        Elimit = self.Elimit * c_list[3]
        #@todo: modify this for the case tension stiffening 
        if e_max <= Epp:
            return 0     
        else:
            return -0.5 * Epp * (- Epp - 2.0*Efp + 2.0*Efp * exp( ((-e_max + Epp) / Efp))) \
                   -0.5 *  e_max * Epp * exp(- ((e_max - Epp) / Efp));

    def get_value(self, e_max, *c_list):
        '''
        Evaluate the integrity of a particular microplane.
        '''
#        print 'x3c_list', c_list
        if len( c_list ) == 0:
            c_list = [1.,1.,1.,1.]
#        print 'x4c_list ', c_list

#        print 'self.Epp used for TensionStiffening:', self.Epp
#        print 'self.Efp used for TensionStiffening:', self.Efp
#        print 'self.Dfp used for TensionStiffening:', self.Dfp
#        print 'self.Elimit used for TensionStiffening:', self.Elimit

        Epp    = self.Epp    * c_list[0]
        Efp    = self.Efp    * c_list[1]
        Dfp    = self.Dfp    * c_list[2]
        Elimit = self.Elimit * c_list[3]
        #
        if e_max <= Epp:
            return 1.0
        elif e_max >= Elimit:
            return 1.e-100
        else:
            return (1-Dfp) * sqrt( Epp / e_max * exp( -(e_max-Epp)/Efp )) + Dfp

    def __call__( self, e_max, *c_list ):
        return self.get_value( e_max, *c_list )

    # Default TraitsUI view
    traits_view = View( Group( Item('Epp'),
                               Item('Efp'),
                               Item('Dfp'),
                               Item('Elimit'),
                               Item( 'mfn', show_label = False),# editor = mfn_editor ),
                               label = 'Damage law',
                               show_border = True
                               ),
                        resizable=True,
                        width=800, height=800)


if __name__ == '__main__':    
    phi_fn = PhiFnQuasiBrittle( Epp = 0.2, Efp = 0.6 )
#    phi_fn = PhiFnTensionStiffening(Epp = 0.2, Efp =0.6, Dfp = 0.2, Elimit = 4.0 )
    phi_fn.configure_traits()