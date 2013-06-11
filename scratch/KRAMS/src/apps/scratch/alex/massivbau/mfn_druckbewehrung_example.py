
from enthought.traits.api import \
    Array, Bool, Callable, Enum, Float, HasTraits, \
    Instance, Int, Trait, Str, Enum, Callable, List, TraitDict, Any, \
    Dict, Property, cached_property, WeakRef, Delegate, \
    ToolbarButton, on_trait_change, Code, Expression, Button
    
from enthought.traits.ui.api import \
    Item, View, HGroup, ListEditor, VGroup, VSplit, Group, HSplit

from enthought.traits.ui.menu import \
    NoButtons, OKButton, CancelButton, Action, CloseAction, Menu, \
    MenuBar, Separator
    #                                 
from mathkit.mfn.mfn_line.mfn_line import MFnLineArray
from mathkit.mfn.mfn_line.mfn_matplotlib_multiline_editor import MFnMatplotlibEditor

#from mathkit.mfn.mfn_line.mfn_chaco_multiline_editor import MFnChacoEditor
#from apps.scratch.faezeh.multiline.mfn_chaco_multiline_editor import MFnChacoEditor
from mfn_druckbewehrung_editor import MFnChacoEditor

from mathkit.mfn.mfn_line.mfn_plot_adapter import MFnPlotAdapter
from numpy import linspace, frompyfunc, vstack, column_stack
from math import sin, cos

class AnalyticalFunction( HasTraits ):
    
    expression1 = Expression('x**2', auto_set = False, enter_set = True )
    expression2 = Expression('x**2.3', auto_set = False, enter_set = True )
    expression3 = Expression('x**2.6', auto_set = False, enter_set = True )
    refresh = Button('redraw')
    def _refresh_fired(self):
        xdata = linspace(0,10,10)
        fneval1 = frompyfunc( lambda x: eval( self.expression1 ), 1, 1 )
        fneval2 = frompyfunc( lambda x: eval( self.expression2 ), 1, 1 )
        fneval3 = frompyfunc( lambda x: eval( self.expression3 ), 1, 1 )
        y1 = fneval1( xdata )
        y2 = fneval2( xdata )
        y3 = fneval3( xdata )
        ydata = column_stack((y1,y2,y3))
        self.mfn.set( xdata = xdata, ydata = ydata )
        self.mfn.data_changed = True
        
    mfn = Instance( MFnLineArray )
    def _mfn_default( self ):
        return MFnLineArray()
    
    @on_trait_change('expression' )
    def update_mfn(self):
        self._refresh_fired()
    
    view_mpl = View( HGroup( Item( 'expression1' ),
                             Item( 'expression2' ),
                             Item( 'expression3' ),
                             Item('refresh' ) ),
                 Item( 'mfn', editor = MFnMatplotlibEditor(), show_label = False ),
                 resizable = True,
                 scrollable = True,
                 height = 0.5, width = 0.5
                    )
    
    view_chaco = View( HGroup( Item( 'expression1' ),
                               Item( 'expression2' ),
                               Item( 'expression3' ), 
                               Item('refresh' ) ),
                 Item( 'mfn', editor = MFnChacoEditor(), resizable = True, show_label = False ),
                 resizable = True,
                 scrollable = True,
                 height = 0.8, width = 0.8
                    )

if __name__ == '__main__':
    fn = AnalyticalFunction()
    fn._refresh_fired()
   # fn.configure_traits( view = "view_mpl", kind = 'nonmodal' )
    fn.configure_traits( view = "view_chaco", kind = 'nonmodal' )
    