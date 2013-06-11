from numpy import array,linspace
from enthought.traits.api import Array, Bool, Callable, Enum, Float, HasTraits, \
                                 Instance, Int, Trait, ToolbarButton, Button, on_trait_change

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
from dacwt import DAC

class MFnWTHandler( Handler ):
    def print_data( self, info ):
        '''
        Eval in a thread.
        '''
        info.object.print_data()

    def save_figure( self, info ):
        '''
        Eval in a thread.
        '''
        info.object.print_data()


class MFnLineArray(DAC):

    # Public Traits
    xdata = Array
    ydata = Array

    def get_value(self,x):
        x2idx = self.xdata.searchsorted(x)
        if x2idx == len(self.xdata):
            x2idx -= 1
        x1idx = x2idx - 1
        x1 = self.xdata[ x1idx ]
        x2 = self.xdata[ x2idx ]
        dx = x2 - x1
        y1 = self.ydata[ x1idx ]
        y2 = self.ydata[ x2idx ]
        dy = y2 - y1
        y = y1 + dy / dx * (x - x1)
        return y

    def get_diff( self, x ):
        x2idx = self.xdata.searchsorted(x)
        if x2idx == len(self.xdata):
            x2idx -= 1
        x1idx = x2idx - 1
        x1 = self.xdata[ x1idx ]
        x2 = self.xdata[ x2idx ]
        dx = x2 - x1
        y1 = self.ydata[ x1idx ]
        y2 = self.ydata[ x2idx ]
        dy = y2 - y1
        return dy / dx
        
    plot_type = Enum("line","scatter")
    dump_button = ToolbarButton('Print data',
                                style = 'toolbar')
    @on_trait_change('dump_button')
    def print_data(self,event = None):
        print repr( self.xdata )
        print repr( self.ydata )

    refresh_button = ToolbarButton('Refresh',
                                   style = 'toolbar')

    @on_trait_change('refresh_button')
    def refresh(self,event = None):
        self.timer_tick()

    # Default TraitsUI view
    view_container = View(Item('container',
                               editor = ComponentEditor(),
                               show_label = False),
                          buttons=NoButtons,
                          handler=MFnWTHandler(),
                          menubar=MenuBar(Menu(Action(name="Print Data",
                                                      action="print_data"),
                                               CloseAction,
                                               name="File")),
                          width=800,
                          height=900,
                          resizable=True)
                      

    traits_view = View( Group( Item("plot_type"),
                               Item('refresh_button', show_label = False ),
                               Item('dump_button', show_label = False),
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
                        handler=MFnWTHandler(),                       
                        menubar=MenuBar(Menu(Action(name="Print data",
                                                    action="print_data"),
                                            Action(name = 'Save figure',
                                                   action='save_figure'),
                                             name="View")),
                        buttons= [OKButton,CancelButton],
                        resizable=True,
                        width=500, height=800)
    
class MFnLineList(MFnLineArray):
    def __init__(self, **kwtraits):
        super(MFnLineList, self).__init__(**kwtraits)
        self.plot_type = 'line'
        self.clear()

    def clear(self):
        self.xdata = array([])
        self.ydata = array([])
        self._xdata = []
        self._ydata = []
        
    def add_pair( self, x, y ):
        self._xdata.append( x )
        self._ydata.append( y )
        
    def set_lists(self, xlist, ylist ):
        self._xdata = xlist
        self._ydata = ylist

    def timer_tick( self, e = None ):
        self.xdata = array( self._xdata )
        self.ydata = array( self._ydata )
        
import wx

class MyApp(wx.PySimpleApp):
    
    def OnInit(self, *args, **kw):

        self.trace = MFnLineList()
        
        self.trace.edit_traits()
        # Set up the timer and start it up
        self.setup_timer()

        from threading import Thread
        eval_thread = Thread( target=self.eval_fn, name='eval')
        eval_thread.start()

        return True

    def setup_timer(self):
        # Create a new WX timer
        timerId = wx.NewId()
        self.timer = wx.Timer(self, timerId)
        
        # Register a callback with the timer event
        self.Bind(wx.EVT_TIMER, self.trace.timer_tick, id=timerId)
        #self.Bind(wx.EVT_TIMER, controller.timer_tick, id=timerId)
        
        return

    def eval_fn( self ):

        # Start up the timer!  We have to tell it how many milliseconds
        # to wait between timer events.  For now we will hardcode it
        # to be 100 ms, so we get 10 points per second.
        self.timer.Start(1000.0, wx.TIMER_CONTINUOUS)

        for i in range( 0, 100 ):
            time.sleep(0.1)
            x = i/10.
            y = x*math.sin(x)
            self.trace.add_pair( x, y )

        self.timer.Stop()
        self.trace.timer_tick()

if __name__ == "__main__":

    #mf = MFnLineArray( xdata = linspace(0,10e-3,10),
    #                    ydata = linspace(0,10e-3,10) )
    #print mf.get_value( 3.5 )
    #mf.configure_traits( )
    #mf.configure_traits( view = view_plot  )
    app = MyApp()
    app.MainLoop()


