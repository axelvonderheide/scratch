from numpy import array,linspace, pi, arange, sin, cos, ones, frompyfunc

from enthought.chaco.api import create_polar_plot

from enthought.enable.example_support import DemoFrame, demo_main

# Enthought library imports
from enthought.enable.api import Window
from enthought.traits.api import false

from enthought.traits.api import Array, Bool, Callable, Enum, Float, HasTraits, \
                                 Instance, Int, Trait, ToolbarButton, Button, on_trait_change, \
                                 Property, cached_property, Range

from enthought.traits.ui.api import Item, View, Group, Handler

from enthought.traits.ui.menu import NoButtons, OKButton, CancelButton, Action, CloseAction, Menu, \
                                     MenuBar, Separator

# Chaco imports
from enthought.chaco.chaco_plot_editor import ChacoPlotItem
# MFnWT imports
from math_func_editor import MFnWTPlotItem

from enthought.enable.component_editor import ComponentEditor

# Chaco imports
from enthought.chaco.example_support import DemoFrame, demo_main, COLOR_PALETTE

from enthought.chaco.api import create_line_plot, add_default_axes, \
                                 add_default_grids, OverlayPlotContainer, \
                                 PlotLabel, VPlotContainer, \
                                 create_scatter_plot, Legend, PlotComponent
from enthought.chaco.tools.api import PanTool, RectZoomTool, SimpleZoom, \
                                       LegendTool, TraitsTool, DragZoom
import time, sys
import math

def radius_fn( theta ):
    if  (theta) < pi/4. or \
        (theta) > 7*pi/4. or \
        ((theta) > 3./4.*pi and (theta) < 5./4.*pi):
        return 0.5
    else:
        return 0.25

vradius_fn = frompyfunc( radius_fn, 1, 1 )

class MFnWTHandler( Handler ):
    def open_data( self, info ):
        info.object.print_data()

    def save_file(self,info):
        sys.exit(0)

    def exit_file(self,info):
        sys.exit(0)

    def init_file(self, window):
        """ Creates a new action. """
        self._window = window
        self.name = "E&xit"

    def perform(self):
        """ Performs the action. """
        self._window.close()

class MFnLineArray(HasTraits):

    numpoints = Int(500000)
    low = Float(0)
    high = Float(2*pi)

    theta  = Array
    def _theta_default(self):
        return arange(self.low, self.high, (self.high-self.low) / self.numpoints)
    
    radius = Property( Array, depends_on = 'theta' )
    @cached_property
    def _get_radius(self):
        return array( vradius_fn( self.theta ), dtype='float_' ) 
    
    plot_type = Enum('polar',"line","scatter")
    
    traits_view = View(MFnWTPlotItem("theta", "radius", 
                                      type_trait="plot_type",
                                      
                                      # Basic axis and label properties
                                      show_label=False,
                                      resizable=True,
                                      orientation="h",
                                      x_label = "Index data",
                                      y_label = "Value data",
                                      
                                      # Plot properties
                                      color = "green",
                                      bgcolor = "white",
                                      
                                      # Specific to scatter plot
                                      marker = "circle",
                                      marker_size = 2,
                                      outline_color = "none",
                                      
                                      # Border, padding properties
                                      border_visible=True,
                                      border_width=1,
                                      padding_bg_color = "lightgray"),
                       buttons= [OKButton,CancelButton],
                       menubar=MenuBar(Menu(Action(name="O&pen..", action="open_data"),
                                            Action(name="S&ave..", action="save_file"),
                                            Action(name="E&xit", action="exit_file"),
                                            name = 'File')),
                       handler = MFnWTHandler,
                       resizable=True,
                       width=500, height=800)


mp = MFnLineArray()
mp.configure_traits()
    
    
    

  