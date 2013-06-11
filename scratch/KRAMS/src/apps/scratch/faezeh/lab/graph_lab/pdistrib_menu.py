from enthought.traits.api import HasTraits, Float, Int, Array, Interface, Tuple, Property, cached_property, \
                                 Any, Enum, Instance, on_trait_change, Dict, false, Button
from enthought.traits.ui.api import Item, View, Group, Handler, HSplit, VSplit, HGroup, HSplit, VGroup
from enthought.traits.ui.menu import NoButtons, OKButton, CancelButton, Action, Menu, MenuBar 
from scipy import stats
from enthought.chaco.example_support import COLOR_PALETTE
from enthought.enable.component_editor import ComponentEditor
# Chaco imports
from enthought.chaco.api import  Plot, ArrayDataSource, create_line_plot, add_default_grids, \
                                 OverlayPlotContainer, PlotLabel, LinearMapper, \
                                 Legend, PlotAxis, add_default_axes, FilledLinePlot, \
                                 DataRange1D, LinePlot, PolygonPlot, DataLabel
from enthought.chaco.tools.api import  PanTool, ZoomTool,  SimpleZoom, \
                                       LegendTool, TraitsTool, BroadcasterTool
# Major library imports
from numpy import  pi, linspace, math, array
import sys

from pdistrib_tool_lab import PDistrib

#from pdistrib_tool_lab import IPDistrib
                                       
class PDistribMenu(HasTraits):

    idistrib = Instance(PDistrib)

    #Exit the form
    def exit_file(self):
        sys.exit(0)

    traits_view = View( 
                        Item('idistrib', style="custom"),

                        menubar=MenuBar(Menu(Action(name ="E&xit",action="exit_file"),
                                             name="View")),
                        buttons= [OKButton,CancelButton],
                        resizable=True,
#                        width=900, height=800
                        )
if __name__ == '__main__':
    idis = PDistrib()
    cl = PDistribMenu(idistrib = idis)
    cl.configure_traits()