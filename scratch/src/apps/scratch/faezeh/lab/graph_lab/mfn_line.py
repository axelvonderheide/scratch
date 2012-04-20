from numpy import array,linspace, trapz, arange
from enthought.traits.api import Array, Bool, Callable, Enum, Float, HasTraits, \
                                 Instance, Int, Trait, ToolbarButton, Button, on_trait_change, \
                                 Property, cached_property

from enthought.traits.ui.api import Item, View, Group, Handler

from enthought.traits.ui.menu import NoButtons, OKButton, CancelButton, Action, CloseAction, Menu, \
                                     MenuBar, Separator

# Chaco imports
from enthought.chaco.chaco_plot_editor import ChacoPlotItem
# MFnWT imports
from mfn_line_editor import MFnWTPlotItem

from enthought.enable.component_editor import ComponentEditor

# Chaco imports
from enthought.chaco.example_support import DemoFrame, demo_main, COLOR_PALETTE

from enthought.chaco.api import create_line_plot, add_default_axes, \
                                 add_default_grids, OverlayPlotContainer, \
                                 PlotLabel, VPlotContainer, \
                                 create_scatter_plot, Legend, PlotComponent
from enthought.chaco.tools.api import PanTool, RectZoomTool, SimpleZoom, \
                                       LegendTool, TraitsTool, DragZoom

import time
import math

class MFnLineArray(HasTraits):

    # Public Traits
    xdata = Array
    def _xdata_default(self):
        '''
        convenience default - when xdata not defined created automatically as
        an array of integers with the same shape as ydata
        '''
        return arange(self.ydata.shape[0])
    
    ydata = Array
    
    plot_type = Enum("line","scatter")
    dump_button = ToolbarButton('Print data',
                                style = 'toolbar')
    @on_trait_change('dump_button')
    def print_data(self,event = None):
        print 'x = ', repr( self.xdata )
        print 'y = ', repr( self.ydata )

    integ_value = Property( Float(), depends_on = 'ydata' ) 
    @cached_property
    def _get_integ_value(self):
        _xdata = self.xdata
        _ydata = self.ydata
        # integral under the stress strain curve
        E_t = trapz(_ydata,_xdata)
        # area of the stored elastic energy  
        U_t = 0.0
        if len(_xdata) != 0:
            U_t = 0.5 * _ydata[-1] * _xdata[-1]        
        return E_t - U_t        

                 

    traits_view = View( Group( Item("plot_type"),
                               Item('dump_button', show_label = False),
                               Item('integ_value', style = 'readonly' ),
                               orientation = 'horizontal' ),
                        MFnWTPlotItem("xdata", "ydata", 
                                      type_trait="plot_type",
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
                        menubar=MenuBar(Menu(Action(name="Print data",
                                                    action="print_data"),
                                            Action(name = 'Save figure',
                                                   action='save_figure'),
                                             name="View")),
                        buttons= [OKButton,CancelButton],
                        resizable=True,
                        width=500, height=800)
    
        

if __name__ == "__main__":

    mf = MFnLineArray( ydata = array([0,1,1.1]) )
  #  print mf.get_value( 5 )
    mf.configure_traits( )


