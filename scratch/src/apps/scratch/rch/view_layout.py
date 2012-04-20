
from enthought.traits.api import \
    HasTraits, Float, Int, List, Array, Interface, String, Tuple, Property, cached_property, \
    Any, Enum, Instance, on_trait_change, Dict, false

from enthought.traits.ui.api import Item, View, Group, Handler, HSplit, VSplit, \
                                    HGroup, HSplit, VGroup, Tabbed

from enthought.traits.ui.menu import NoButtons, OKButton, CancelButton, Action, CloseAction, Menu, \
                                     MenuBar, Separator

from scipy import stats
from numpy import linspace
import math

# Plotting 

from enthought.chaco.example_support import COLOR_PALETTE

from enthought.enable.component_editor import \
    ComponentEditor
from enthought.chaco.api import \
    Plot, AbstractPlotData, ArrayPlotData, \
    ArrayDataSource
from enthought.chaco.tools.api import \
    PanTool, ZoomTool
    
from math import sqrt

# Major library imports
from numpy import arange, fabs, pi, sin, array, sqrt
from scipy.special import jn

# Chaco imports
from enthought.chaco.api import \
    create_line_plot, add_default_axes, add_default_grids, \
    OverlayPlotContainer, PlotLabel, VPlotContainer, \
    create_scatter_plot, Legend, PlotComponent, PlotAxis, \
    PlotAxis, FilledLinePlot, \
    DataRange1D, LinePlot, BarPlot, PolygonPlot, LinearMapper

from enthought.chaco.tools.api import PanTool, RectZoomTool, SimpleZoom, \
                                       LegendTool, TraitsTool, BroadcasterTool
                          

from scipy.optimize import fsolve

pdist = {}

pdist['norm'] = stats.norm
pdist['weibull'] = stats.weibull_min
pdist['uniform'] = stats.uniform
pdist['gamma'] = stats.gamma

                          
class PDistrib(HasTraits):

    #------------------------------------------------------------------------
    # Distribution parameters
    #------------------------------------------------------------------------
    loc   = Float( 1.0, auto_set = False, enter_set = True )
    scale = Float( 1.0, auto_set = False, enter_set = True )
    shape = Float( 1.0, auto_set = False, enter_set = True )

    distr_type = Enum('norm','gamma','weibull','uniform', 'cos')

    mean = Float(1.0)
    variance = Float(2.0)
    kurtosis = Float(0, auto_set = False, enter_set = True )

    traits_view = View( HSplit( Group( 
                                           Item('distr_type'),                                       
                                           Item('scale'),
                                           Item('shape'),                                           Item('loc'),
                                           Item('mean'),
                                           Item('variance'),
                                           Item('skew'),
                                           label = 'parameters',
                                           scrollable = True
                                            ),
                                        Group(
                                           Item('quantile'),
                                           Item('n_segments'),
                                           label = 'plot ctrls',
                                           scrollable = True
                                            ),
                                   scrollable = True,
                                    id = 'stats.pdistrib.hsplitter'
                                ),
                        menubar=MenuBar(Menu(Action(name="Print data",
                                                        action="print_data"),
                                                Action(name = 'Save figure',
                                                   action='save_figure'),
                                             name="View")),
                        dock = 'tab',
                        #id = 'stats.pdistrib.pdistrib2',                                       
                        buttons= [OKButton,CancelButton],
                        resizable=True,
                        width=0.4, height=0.4
                        )


if __name__ == '__main__':
    normal = PDistrib( distr_type = 'uniform' )
    normal.configure_traits()
