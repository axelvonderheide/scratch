
from enthought.traits.api import \
    Array, Bool, Callable, Enum, Float, HasTraits, \
    Instance, Int, Trait, Str, Enum, Callable, List, TraitDict, Any, \
    Dict, Property, cached_property, WeakRef, Delegate, \
    ToolbarButton, on_trait_change

from enthought.traits.ui.api import \
    Item, View, HGroup, ListEditor, VGroup, VSplit, Group, HSplit, Spring

from enthought.traits.ui.menu import \
    NoButtons, OKButton, CancelButton, Action, CloseAction, Menu, \
    MenuBar, Separator

from numpy import float_, zeros, arange, array, copy, frompyfunc, vstack
import pprint
from ibvpy.api import RTrace

from mathkit.mfn.mfn_line.mfn_line import MFnLineArray
from mathkit.mfn.mfn_line.mfn_matplotlib_editor import MFnMatplotlibEditor
from mathkit.mfn.mfn_line.mfn_plot_adapter import MFnPlotAdapter

class RTraceGraph( RTrace ):
    '''
    Collects two response evaluators to make a response graph.
    
    The supplied strings for var_x and var_y are used to locate the rte in 
    the current response manager. The bind method is used to navigate to 
    the rte and is stored in here as var_x_eval and var_y_val as Callable 
    object.
    
    The request for new response evaluation is launched by the time loop 
    and directed futher by the response manager. This method is used solely 
    for collecting the data, not for their visualization in the viewer.
    
    The timer_tick method is invoked when the visualization of the Graph
    should be synchronized with the actual contents.
    '''

    label = Str( 'RTraceGraph' )
    var_x = Str( '' )
    var_x_eval = Callable
    idx_x_arr = Array
    idx_x = Int( -1, enter_set = True, auto_set = False )
    var_y = Str( '' )
    var_y_eval = Callable
    idx_y_arr = Array
    idx_y = Int( -1, enter_set = True, auto_set = False )
    transform_x = Str( enter_set = True, auto_set = False )
    transform_y = Str( enter_set = True, auto_set = False )

    trace = Instance( MFnLineArray )
    def _trace_default( self ):
        return MFnLineArray()

    print_button = ToolbarButton( 'Print Values',
                                   style = 'toolbar' )

    @on_trait_change( 'print_button' )
    def print_values( self, event = None ):
        print 'x:\t', self.trace.xdata, '\ny:\t', self.trace.ydata

    view = View( VSplit( VGroup( HGroup ( VGroup( HGroup( Spring(), Item( 'var_x', style = 'readonly' ),
                                                         Item( 'idx_x', show_label = False ) ),
                                                  Item( 'transform_x' ) ),
                                          VGroup( HGroup( Spring(), Item( 'var_y', style = 'readonly' ),
                                                         Item( 'idx_y', show_label = False ) ),
                                                         Item( 'transform_y' ) ),
                                          VGroup( 'update_on', 'clear_on' ) ),
                                 HGroup( Item( 'refresh_button', show_label = False ),
                                         Item( 'print_button', show_label = False ) ),
                                ),
                         Item( 'trace@', \
                              editor = MFnMatplotlibEditor( \
                                        adapter = MFnPlotAdapter( var_x = 'var_x',
                                                                  var_y = 'var_y',
                                                                  min_size = ( 100, 100 ),
                                                                  max_size = ( 350, 350 ) ) ),
                              show_label = False, resizable = True ),
                         ),
                 buttons = [OKButton, CancelButton],
                 resizable = True,
                 scrollable = True,
                 height = 0.5, width = 0.5 )

    _xdata = List( Array( float ) )
    _ydata = List( Array( float ) )

    def bind( self ):
        '''
        Locate the evaluators
        '''
        self.var_x_eval = self.rmgr.rte_dict.get( self.var_x, None )
        if self.var_x_eval == None:
            raise KeyError, 'Variable %s not present in the dictionary:\n%s' % \
                            ( self.var_x, self.rmgr.rte_dict.keys() )

        self.var_y_eval = self.rmgr.rte_dict.get( self.var_y, None )
        if self.var_y_eval == None:
            raise KeyError, 'Variable %s not present in the dictionary:\n%s' % \
                            ( self.var_y, self.rmgr.rte_dict.keys() )


    def setup( self ):
        self.clear()

    def add_current_values( self, sctx, U_k, *args, **kw ):
        '''
        Invoke the evaluators in the current context for the specified control vector U_k.
        '''
        x = self.var_x_eval( sctx, U_k, *args, **kw )
        y = self.var_y_eval( sctx, U_k, *args, **kw )

        self.add_pair( x.flatten(), y.flatten() )

    def add_pair( self, x, y ):
        self._xdata.append( copy( x ) )
        self._ydata.append( copy( y ) )

    @on_trait_change( 'idx_x,idx_y' )
    def redraw( self, e = None ):
        if ( ( self.idx_x < 0 and len( self.idx_x_arr ) == 0 ) or
             ( self.idx_y < 0 and len( self.idx_y_arr ) == 0 ) or
             self._xdata == [] or
             self._ydata == [] ):
            return
        #
        if len( self.idx_x_arr ) > 0:
            print 'x: summation for', self.idx_x_arr
            xarray = array( self._xdata )[:, self.idx_x_arr].sum( 1 )
        else:
            xarray = array( self._xdata )[:, self.idx_x]

        if len( self.idx_y_arr ) > 0:
            print 'y: summation for', self.idx_y_arr
            yarray = array( self._ydata )[:, self.idx_y_arr].sum( 1 )
        else:
            yarray = array( self._ydata )[:, self.idx_y]

        if self.transform_x:
            def transform_x_fn( x ):
                '''makes a callable function out of the Str-attribute
                "transform_x". The vectorised version of this function is 
                then used to transform the values in "xarray". Note that
                the function defined in "transform_x" must be defined in
                terms of a lower case variable "x".
                '''
                return eval( self.transform_x )
            xarray = frompyfunc( transform_x_fn, 1, 1 )( xarray )

        if self.transform_y:
            def transform_y_fn( y ):
                '''makes a callable function out of the Str-attribute
                "transform_y". The vectorised version of this function is 
                then used to transform the values in "yarray". Note that
                the function defined in "transform_y" must be defined in
                terms of a lower case variable "y".
                '''
                return eval( self.transform_y )
            yarray = frompyfunc( transform_y_fn, 1, 1 )( yarray )

        self.trace.xdata = array( xarray )
        self.trace.ydata = array( yarray )
        self.trace.data_changed = True

    def timer_tick( self, e = None ):
        # @todo: unify with redraw
        pass

    def clear( self ):
        self._xdata = []
        self._ydata = []
        self.trace.clear()
        self.redraw()

class RTraceArraySnapshot( RTrace ):
    '''
    Plot the current value of the array along the x_axis
    
    Used currently for plotting the integrity factor over microplanes 
    '''
    var = Str( '' )
    var_eval = Callable
    idx = Int( -1 )

    trace = Instance( MFnLineArray )
    def _trace_default( self ):
        return MFnLineArray()

    view = View( HSplit( VGroup ( VGroup( 'var' ),
                                  VGroup( 'update_on', 'clear_on' ) ),
                         Item( 'trace@', style = 'custom',
                              editor = MFnMatplotlibEditor( \
                                        adapter = MFnPlotAdapter( var_y = 'var' ) ),
                              show_label = False, resizable = True ),
                         ),
                 buttons = [OKButton, CancelButton],
                 resizable = True )

    def bind( self ):
        '''
        Locate the evaluators
        '''
        self.var_eval = self.rmgr.rte_dict[self.var]

    def setup( self ):
        pass

    def add_current_values( self, sctx, U_k, *args, **kw ):
        '''
        Invoke the evaluators in the current context for the specified control vector U_k.
        '''
        self.y = array( self.var_eval( sctx, U_k, *args, **kw ), dtype = 'float' )
        self.x = arange( 0, len( self.y ), dtype = 'float' )
        self.redraw()

    def timer_tick( self, e = None ):
        self.redraw()

    def redraw( self ):
        self.trace.xdata = self.x
        self.trace.ydata = self.y

    def clear( self ):
        # @todo: 
        self.x = zeros( ( 0, ), dtype = 'float_' )
        self.y = zeros( ( 0, ), dtype = 'float_' )

if __name__ == '__main__':

    rm1 = RTraceGraph( name = 'rte 1',
                       idx_x = 0,
                       idx_y = 0,
                       transform_x = '-x' )
#                       transform_x = lambda x: -x )

    rm1.add_pair( array( [0.0, 0.0] ), array( [0.0, 0.0] ) )
    rm1.add_pair( array( [1.0, 1.0] ), array( [1.0, 1.5] ) )
    rm1.add_pair( array( [2.0, 1.5] ), array( [2.0, 2.0] ) )
    rm1.redraw()
    rm1.configure_traits()
