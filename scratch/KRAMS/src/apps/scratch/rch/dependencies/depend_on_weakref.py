
# Traits imports
from enthought.traits.api import \
     Array, Bool, Callable, Enum, Float, HasTraits, \
     Instance, Int, Trait, Range, HasTraits, on_trait_change, Event, \
     implements, Dict, Property, cached_property, Delegate, WeakRef, String, \
     Constant, List

# Traits UI imports
from enthought.traits.ui.api import \
     Item, View, HSplit, VSplit, VGroup, Group, Spring, TabularEditor
from enthought.traits.ui.menu \
    import OKButton, CancelButton
from enthought.traits.ui.tabular_adapter \
    import TabularAdapter


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

# Numpy imports
from numpy import \
     array, ones, zeros, outer, inner, transpose, dot, frompyfunc, \
     fabs, sqrt, linspace, vdot, identity, tensordot, \
     sin as nsin, meshgrid, float_, ix_, \
     vstack, hstack, sqrt as arr_sqrt

from math import pi as Pi, cos, sin, exp, sqrt as scalar_sqrt

from scipy.linalg import eig, inv

from core.tstepper import \
     TStepper as TS

from core.mats_eval import IMATSEval, MATSEval

from core.rv import RTrace, RTraceGraph, RTraceArraySnapshot

from subsid.mathkit import MFnNDGrid

#----------------------------------------------------------------------------------
#                                     VariedParam
#----------------------------------------------------------------------------------
class VariedParam(HasTraits):
    """
    Association between the spatial function and material parameter.
    """
    mats_eval = WeakRef( IMATSEval )

    varname = String
    
    reference_value = Property( Float, depends_on = 'mats_eval' )
    @cached_property
    def _get_reference_value(self):
        return getattr( self.mats_eval, self.varname )
    
    switch = Enum('constant','varied')
    
    spatial_fn = Instance( MFnNDGrid )
    
    variable = Property( Float )
    def _get_variable(self):
        return getattr( self.mats_eval, self.varname )
    def _set_variable(self, value):
        setattr( self.mats_eval, self.varname, value )
        
    def adjust_for_sctx( self, sctx ):
        if self.switch == 'varied':
            X_pnt = fets_eval.get_X_pnt( sctx )
            coeff = self.spatial_fn( X_pnt )
            self.variable = self.reference_value * coeff
        
    traits_view = View( Group( Item( 'varname', style = 'readonly', show_label = False ),
                               Item( 'reference_value', style = 'readonly'),
                               Item( 'switch', style = 'custom', show_label = False ),
                               Item( 'spatial_fn', style = 'custom', show_label=False, resizable = True ) ),
                        resizable = True,
                        height=800 )

#----------------------------------------------------------------------------------
# Tabular Adapter Definition 
#----------------------------------------------------------------------------------
class VariedParamAdapter ( TabularAdapter ):

    columns = [ ( 'Name',         'varname' ), 
                ( 'Variable',     'variable' ) ]
                
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

#---------------------------------------------------------------------------
# Material Proxy to include variable parameters
#---------------------------------------------------------------------------
class MATSProxy( MATSEval ):
    '''
    Material model with spatially varying material parameters.
    
    @author: rch    
    The material model works as a proxy delegating the standard 
    functionality to the associated material model @param mats_eval.
    
    In the first step, the proxy identifies the material parameters 
    of the mats_eval. The identification is performed by scanning the 
    @param mats_eval for float traits using the method identify parameters.
    The scanning happens on demand when the @param 
    varpars list is accessed.

    For each variable an instance of the VarPar object gets constructed 
    to handle the spatial variation. By default, the VarPar has the 
    switch set to constant so that the proxy has no effect. 
    
    Spatial variation of a parameter within the VarPar instance 
    can be activated both in the script during the model construction 
    and interactively thorough the user interface. Both these ways are
    now briefly described to illuminate the internal workings 
    of the proxy.
    
    Scripting interface
    ===================
    There is a quite number of depdencies involved in associating a
    spatial profile to a material variable. Since this implementation
    is meant generic to be applicable to any material model, 
    these dependcies must be captured in a transparent way. Let us
    identify the dependencies using a particular example.
    
    *** initial state *** 

    proxy has an associated default material model. 

    (@note: actually, it would be reasonable to reflect the usage 
    context of the proxy, i.e. reuse the mats_eval used within of the element 
    element using the mats. Probably the spatial variation should intervene
    at a global level and hook up the spatial context anywhere upon the 
    change of the spatial coordinate. In other words, 
    the adjustment of a parameter in a particular time stepper would be done
    whenever the spatial coordinate gets changed.  
 
    However, handling of these association 
    is not the current priority. It can be solved by subclassing 
    the proxy for the five supported dimensions, i.e. 1d, 1d5, 2d, 2d5 and 3d)
     
    upon access to varpar: identify the parameters of 
    @param mats_eval and construct the VarPar instances. They are all
    set to constant.
    
    *** assign a spatial profile ***
      
    The variables must be accessible via their keywords. This means they are 
    manged within a dictionary. For an instance of the form you can issue
    
    mp = MATSProxy( mats_eval = MATS1DElastic )
    mp.varpars]'E'].spatial_fn.x_mins = [0,0,0]
    mp.varpars]'E'].spatial_fn.x_mins = [0,0,0]
    mp.varpars]'E'].spatial_fn.x_mins = [0,0,0]

    mp['E'].spatial_fn.shape = 500
    mp['E'].spatial_fn.set_values_in_box( 10, [-1,-1,-1], [8,8,8] )
    
    mp['E'].spatial_fn.shape = 800
    mp['E'].spatial_fn.shape = 800
    
    @todo The spatial variation can be handled at the level 
    of the tstepper     Then, the identification of parameters 
    can be performed for all integration levels. Further, 
    the spatial binding can be taken into account.
    i.e., the spatial profile is specified for an existing 
    domain with particular bounding box. The VaPars are registered 
    in the spatial context and the adjustment 
    of the material parameters is invoked dynamically, 
    whenever the @param r_pnt gets change. That means the correct 
    parameters are available both for 
    the iteration and for response tracing.
    
    The extensions needed - 
    
    1) the recursive list of sub time steppers during the
        parameter identification in tstepper_eval
    
    2) spatial binding to the geometric model as a 
        separate step (identification of the bounding box)
        
    3) spatial context is launched before starting the computation.
        at this point, the varpar_dict must be registered within
        the spatial context.  
        
    3) spatial context must be hooked with the call to 
       the adjust state variables. 
       
    4) any change in the tstepper structure must 
       cause reconstruction of the varpar list. Only those 
       varpars should be modified that are really affected.
       (the need for change notification mechanism within the 
       tstepper hierarchy)
    '''
    implements( IMATSEval )

    mats_eval_type = Enum('MATS1DElastic',
                          'MATS1DDamage',
                          'MA2DCompositeMicroplaneDamage')
    mats_eval = Property( Instance( IMATSEval ), depends_on = 'mats_eval_type' )
    @cached_property
    def _get_mats_eval(self):
        return eval( self.mats_eval_type + '()' )

    #-----------------------------------------------------------------------------------------------
    # Management of spatially varying parameters depending on the value of mats_eval
    #-----------------------------------------------------------------------------------------------
    varpars = Property( List, depends_on = 'mats_eval' )
    @cached_property
    def _get_varpars( self ):
        '''
        reset the varpar list according to the current mats_eval object.
        '''
        params = self.mats_eval.identify_parameters()
        varset = []
        for key, par in params.items():
            par_val = getattr( self.mats_eval , key )
            varset.append (
                    VariedParam(mats_eval = self.mats_eval,
                                varname = key ) )
#                                spatial_fn = 1 ) )
        return varset

    @on_trait_change('varpars')
    def set_current_varpar(self):
        self.current_varpar = self.varpars[0]
        
    # variable selectable in the table of varied params (just for viewing) 
    current_varpar = Instance( VariedParam )
    def _current_varpar_default(self):
        return self.varpars[0]
        
    #---------------------------------------------------------------------------------------------
    # View specification
    #---------------------------------------------------------------------------------------------
    traits_view = View( Group(
                              Item('mats_eval_type', show_label = False, style = 'custom' ),
                              Item('mats_eval', show_label = False, style = 'custom' ),
                              label = 'Material_model'
                              ),
                        Group(
                              HSplit( 
                                    Item('varpars', show_label = False, editor = varpar_editor ),
                                    Item('current_varpar', show_label = False, style = 'custom', resizable = True),
                                     ),
                              label = 'Spatially varied parameters'
                              ),
                        width = 0.8,
                        height = 0.8,
                        resizable = True )


#--------------------------------------------------------------------------------
# Example 
#--------------------------------------------------------------------------------

from core.tloop import TLoop, TLine
from core.bc import BCDof
from core.ibvp_solve import IBVPSolve as IS
from lib.mats.mats1D.mats1D_elastic.mats1D_elastic import MATS1DElastic  
from lib.mats.mats1D.mats1D_damage.mats_damage1d import MATS1DDamage  
from lib.mats.mats2D.mats2D_cmdm.mats_mdm2d import MA2DCompositeMicroplaneDamage
                                       
if __name__ == '__main__':
    # tseval for a material model
    #
    tseval  = MATSProxy(  mats_eval_type = 'MATS1DElastic' )
    tseval.varpars
    tseval.configure_traits()
    