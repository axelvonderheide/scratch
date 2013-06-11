
from numpy import array,linspace, pi, arange, sin, cos, ones, frompyfunc, where, hstack

from enthought.chaco.api import create_polar_plot

from enthought.enable.example_support import DemoFrame, demo_main

# Enthought library imports
from enthought.enable.api import Window
from enthought.traits.api import false
from enthought.traits.api import Delegate
from enthought.traits.api import Array, Bool, Callable, Enum, Float, HasTraits, \
                                 Instance, Int, Trait, ToolbarButton, Button, on_trait_change, \
                                 Property, cached_property, Range, Instance, Str

from enthought.traits.ui.api import Item, View, Group, Handler

from enthought.traits.ui.menu import NoButtons, OKButton, CancelButton, Action, CloseAction, Menu, \
                                     MenuBar, Separator

from enthought.enable.component_editor import ComponentEditor

from enthought.chaco.example_support import DemoFrame, demo_main, COLOR_PALETTE

from enthought.chaco.api import HPlotContainer, create_line_plot, add_default_axes, \
                                 add_default_grids, OverlayPlotContainer, \
                                 PlotLabel, VPlotContainer, \
                                 create_scatter_plot, Legend, PlotComponent
from enthought.chaco.tools.api import PanTool, RectZoomTool, SimpleZoom, \
                                       LegendTool, TraitsTool, DragZoom

from mfn_polar_editor import MFnPolarPlotItem, MFnPolarPlotEditor

import enthought.units as units
from enthought.units.angle import degree, radian

import time, sys
import math


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

class MFnPolar(HasTraits):

    numpoints = Int(80)
    low = Float(0.)
    high = Float(2*pi)
    # arguments to configure the Polar function to be plotted
    alpha                = Range( 0., pi  , 0., auto_set = False)
    delta_alpha          = Range( 0., pi/2, 0., auto_set = False)
    delta_trans          = Range( 0., pi  , 0., auto_set = False)
    phi_residual         = Float( 0.65)
    phi_quasibrittle     = Float( -0.25)
    strech_residual      = Range( 0., 1 , 0., auto_set = False)
    strech_quasibrittle  = Range( 0., 1 , 0., auto_set = False)
#    # arguments to configure the plot
#    plotrange_min        = Float( 0.)
#    plotrange_max        = Float( 1.)
#    frac_noplot          = Range( 0., 1 , 0.3, auto_set = False)
#    # @todo: put this into the status bar
#    info                 = Str('')


    theta  = Array
    def _theta_default(self):
        theta_arr = arange(self.low, self.high, (self.high-self.low) / self.numpoints) 
        # add theta=0 to the end of theta-array
        theta_zero = array([0])
        return hstack([theta_arr,theta_zero])
         
    
    # function to be plotted (returns radius_value for a given theta_value)
    def radius_fn(self, theta, alpha, delta_alpha, delta_trans, strech_residual, strech_quasibrittle, phi_residual, phi_quasibrittle):
        #1st quadrant
        if ((theta-alpha) >= 0. and (theta-alpha) <= (pi/2)):
           theta_tilde = theta-alpha 
        #2nd quadrant   
        elif ((theta-alpha) <= pi and (theta-alpha) >= pi/2):
             theta_tilde = pi-(theta-alpha) 
        #3rd quadrant positive
        elif ((theta-alpha) >= pi and (theta-alpha) <= 3*pi/2):
            theta_tilde = theta-alpha-pi
        #3rd quadrant negative  
        elif ((theta-alpha) >= -pi and (theta-alpha) <= -pi/2):
            theta_tilde = theta-alpha+pi
        #4th quadrant positive
        elif ((theta-alpha) <= 2*pi and (theta-alpha) >= 3*pi/2):
            theta_tilde = (2*pi)-(theta-alpha)
        #4th quadrant negative        
        elif ((theta-alpha) <= 0. and (theta-alpha) >= -pi/2):
            theta_tilde = -(theta-alpha)
    
        ### Definition of function to be plotted in the range of 0 and Pi/2:
        _phi_residual       =  phi_residual + (1-phi_residual) * strech_residual
        _phi_quasibrittle   =  phi_quasibrittle + (1-phi_quasibrittle) * strech_quasibrittle
       
        # constant values with linear transition function:
        # (for delta_alpha = 0 the transition function is evaluated)
        if abs(theta_tilde) < delta_alpha:
            radius_fn = _phi_residual 
        elif abs(theta_tilde) >= delta_alpha and abs(theta_tilde) < delta_alpha+delta_trans:
            radius_fn = (_phi_residual - \
                         ((theta_tilde-delta_alpha)*(_phi_residual - _phi_quasibrittle)/(delta_trans)))
        else:
            radius_fn = _phi_quasibrittle
    
        return radius_fn
 
    
    radius = Property( Array, depends_on = 'theta,alpha,delta_alpha,delta_trans,frac_noplot,strech_residual,strech_quasibrittle,phi_residual,phi_quasibrittle ' )
    @cached_property
    def _get_radius(self):
        vradius_fn = frompyfunc( self.radius_fn, 8, 1 )
        return array( vradius_fn( self.theta, self.alpha, self.delta_alpha, self.delta_trans, self.strech_residual, \
                                  self.strech_quasibrittle, self.phi_residual, self.phi_quasibrittle ), dtype='float_' )


    def __call__(self, theta_value):
        # return a single value for the specified theta_value 
        radius_value = self.radius_fn(theta_value, self.alpha, self.delta_alpha, self.delta_trans, self.strech_residual, self.strech_quasibrittle, self.phi_residual, self.phi_quasibrittle)
        return radius_value

    plot_type = Enum('polar')
    
    
    radiusplot = MFnPolarPlotItem("theta", ["radius"],
                                  type_trait="plot_type",
                                  
                                  # Basic axis and label properties
                                  show_label=False,
                                  resizable=True,
                                  orientation="h",
                                  x_label = "Index data",
                                  y_label = "Value data",
                                  
                                  # Plot properties
                                  color = "green",
                                  bgcolor = "lightyellow",
                                  
                                  # Specific to scatter plot
                                  marker = "circle",
                                  marker_size = 2,
                                  outline_color = "none",
                                  
                                  # Border, padding properties
                                  border_visible=True,
                                  border_width=1,
                                  padding_bg_color = "lightgray")

    
    radius_min = Property( Float, depends_on = 'current_theta,alpha,delta_alpha,delta_trans,strech_residual,strech_quasibrittle' )
    @cached_property
    def _get_radius_min(self):
        r_min = self.radius.min()
        return r_min

    radius_max = Property( Float, depends_on = 'current_theta,alpha,delta_alpha,delta_trans,strech_residual,strech_quasibrittle' )
    @cached_property
    def _get_radius_max(self):
        r_max = self.radius.max()
        return r_max




    current_theta = Float( 0.0 )

    current_radius = Property( Float, depends_on = 'current_theta,alpha,delta_alpha,delta_trans,strech_residual,strech_quasibrittle' )
    @cached_property
    def _get_current_radius(self):
        return self.__call__( self.current_theta )

       
    traits_view = View(Item("alpha"),
                       Item("delta_alpha"),
                       Item("delta_trans"),
                       Item("phi_residual"),
                       Item("phi_quasibrittle"),
                       Item("strech_residual"),
                       Item("strech_quasibrittle"),
#                       Item("plotrange_min"),
#                       Item("plotrange_max"),
#                       Item("frac_noplot"),
                       radiusplot,
                       radiusplot2,
                       Item('current_theta'),
                       Item('current_radius', style = 'readonly' ),
                       Item('radius_min', style = 'readonly' ),
                       Item('radius_max', style = 'readonly' ),
#                       Item('info', style = 'readonly' ),
                       buttons= [OKButton,CancelButton],
                       menubar=MenuBar(Menu(Action(name="O&pen..", action="open_data"),
                                            Action(name="S&ave..", action="save_file"),
                                            Action(name="E&xit", action="exit_file"),
                                            name = 'File')),
                    
                       handler = MFnWTHandler,
                       resizable=True,
                       width=700, height=800)
       
if __name__ == '__main__':
#    mp = MFnPolar( alpha = 0.1, delta_alpha = 0.2, delta_trans = 0.3 )
#    mp(0.224)
#    print 'mp(0.224)',mp(0.224) 
    mp = MFnPolar()
    mp.configure_traits()

