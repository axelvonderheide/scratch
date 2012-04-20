
# Imports:
from enthought.traits.api \
    import HasTraits, HasPrivateTraits, Array, Tuple, Float, Enum, on_trait_change
    
from enthought.traits.ui.api \
    import View, HSplit, VSplit, Item, TupleEditor, Group, Handler
    
from enthought.enable.component_editor import ComponentEditor

from enthought.traits.ui.menu import NoButtons, OKButton, CancelButton, Action, CloseAction, Menu, \
                                     MenuBar
                                     
from mfn_line_editor import MFnWTPlotItem

from numpy import linspace, array

from math import cos

class X2COSX(HasTraits):
    ''' Simple class emulating a function x**2 * cos(x)
    '''

    a = Float(1.)
    b = Float(2.)
    c = Float(3.)


    xdata = Array
    ydata = Array

    @on_trait_change('a,b,c')
    def refresh(self):
        self.xdata = linspace(0.,10.,100)
        self.ydata = array([ self(x) for x in self.xdata ])
        print 'refreshed'
    
    traits_view = View( Group( Item( name = 'a'), 
                               Item( name = 'b'),
                               Item( name = 'c' ),
                               orientation = 'horizontal' ),
                        MFnWTPlotItem("xdata", "ydata", 
                                      type="line",
                                      # Basic axis and label properties
                                      show_label=False,
                                      resizable=True,
                                      orientation="h",
                                      x_label = "X",
                                      y_label = "Y",
                                      title = '',
                                      
                                      # Plot properties
                                      color = "blue",
                                      bgcolor = "white",
                                      
                                      # Specific to scatter plot
                                      marker = "circle",
                                      marker_size = 5,
                                      outline_color = "none",
                                      
                                      # Border, padding properties
                                      border_visible=True,
                                      border_width=0,
                                      padding_bg_color = "white"),
                        buttons= [OKButton,CancelButton],
                        resizable=True,
                        width=500, height=800)
    
    def __call__(self, x):
        return self.a*x**2 * cos(self.b*x ) ** self.c
    
    
if __name__ == '__main__':
    fx = X2COSX()
    print fx(2.4)
    fx.refresh()
    fx.configure_traits()
    
