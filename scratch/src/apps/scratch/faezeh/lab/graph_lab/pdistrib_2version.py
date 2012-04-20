
from enthought.traits.api import \
    HasTraits, Float, Int, List, Array, Interface, String, Tuple, Property, cached_property, \
    Any, Enum, Instance, on_trait_change, Dict, false

from enthought.traits.ui.api import Item, View, Group, Handler, HSplit, VSplit, HGroup, HSplit, VGroup

from enthought.traits.ui.menu import NoButtons, OKButton, CancelButton, Action, CloseAction, Menu, \
                                     MenuBar, Separator

from scipy import stats

from enthought.chaco.example_support import COLOR_PALETTE

from enthought.enable.component_editor import ComponentEditor

# Chaco imports
from enthought.chaco.api import     Plot, AbstractPlotData, ArrayPlotData, ArrayDataSource, \
                                    create_line_plot, add_default_axes, add_default_grids, \
                                    OverlayPlotContainer, PlotLabel, VPlotContainer, LinearMapper, \
                                    create_scatter_plot, Legend, PlotComponent, PlotAxis, FilledLinePlot, \
                                    DataRange1D, LinePlot, BarPlot, TextBoxOverlay, DataLabel, LabelAxis

from enthought.chaco.tools.api import  PanTool, ZoomTool, RectZoomTool, SimpleZoom, \
                                       LegendTool, TraitsTool, BroadcasterTool, DataLabelTool

# Major library imports
from numpy import arange, fabs, pi, sin, linspace, math, array
from scipy.special import jn


 
                                       
class IPDistrib(Interface):
    pass

class PDistrib(HasTraits):

    implements   = IPDistrib

    #------------------------------------------------------------------------
    # Distribution parameters
    #------------------------------------------------------------------------
    loc   = Float( 1.0, auto_set = False, enter_set = True )
    scale = Float( 1.0, auto_set = False, enter_set = True )
    shape = Float( 1.0, auto_set = False, enter_set = True )

    distr_type = Enum('norm','gamma','weibull_min','uniform')

    # Actual distribution implementation (scipy object)
    distr = Property( Any, depends_on = 'loc,scale,shape,distr_type' )
    @cached_property
    def _get_distr(self):
        return stats.__dict__[self.distr_type]( self.shape, loc = self.loc, scale = self.scale )

    #------------------------------------------------------------------------
    # Methods calculating the statistical moments
    #------------------------------------------------------------------------
    stats = Property( Tuple( Float ), depends_on = 'loc,scale,shape,distr_type' )
    @cached_property
    def _get_stats(self):
        print 'refreshing stats'
        return self.distr.stats()

    mean = Property( Float, depends_on = 'loc,scale,shape,distr_type' )
    @cached_property
    def _get_mean(self):
        return self.stats[0]

    
    standard_deviation = Property( Float, depends_on = 'loc,scale,shape,distr_type' )
    @cached_property
    def _get_standard_deviation(self):
        return math.sqrt( self.stats[1] )

    skew = Property( Float, depends_on = 'loc,scale,shape,distr_type' )
    @cached_property
    def _get_skew(self):
        if len( self.stats ) > 2:
            return self.stats[2]
        else:
            return 0

    kurtosis = Property( Float, depends_on = 'loc,scale,shape,distr_type' )
    @cached_property
    def _get_kurtosis(self):
        if len( self.stats ) > 3:
            return self.stats[3]
        else:
            return 0

    #------------------------------------------------------------------------
    # Methods preparing visualization
    #------------------------------------------------------------------------
    quantile = Float(0.00001)
    range = Property( Tuple( Float ), depends_on = 'loc,scale,shape,distr_type,quantile' )
    @cached_property
    def _get_range(self):
        return (self.distr.ppf( self.quantile ), self.distr.ppf( 1-self.quantile ) )

    n_segments = Int(100)

    x_array = Property( Array( Float ), 
                        depends_on = 'loc,scale,shape,distr_type,quantile,n_segments'  )
    @cached_property
    def _get_x_array(self):
        return linspace( self.range[0], self.range[1], self.n_segments )
    
    pdf_array = Property( Array( Float ), 
                          depends_on = 'loc,scale,shape,distr_type,quantile,n_segments' )
    @cached_property
    def _get_pdf_array(self):
        return self.distr.pdf( self.x_array )

    cdf_array = Property( Array( Float ), 
                          depends_on = 'loc,scale,shape,distr_type,quantile,n_segments' )
    @cached_property
    def _get_cdf_array(self):
        return self.distr.cdf( self.x_array )

    plot_container = Instance( OverlayPlotContainer )
    def _plot_container_default(self):
        container = OverlayPlotContainer( padding = 60, fill_padding = False,
                                  bgcolor = "white", use_backbuffer=True)
        # Plot some distribution functions
        plots = {}
        broadcaster = BroadcasterTool()
    
    #""" Plot 
        
#        view = DataView(border_visible = True)
#
        index = ArrayDataSource(self.x_array)
        value = ArrayDataSource(self.pdf_array, sort_order="none")

        index_range = DataRange1D()
        index_range.add(index)
        index_mapper = LinearMapper(range=index_range)
    
        value_range = DataRange1D()
        value_range.add(value)
        value_mapper = LinearMapper(range=value_range)

        pdf_plot = FilledLinePlot(index = index, value = value,
                                    index_mapper = index_mapper,
                                    value_mapper = value_mapper,
                                    edge_color = tuple(COLOR_PALETTE[0]),
                                    face_color = "paleturquoise",
                                    bgcolor = "white",
                                    border_visible = True)

        add_default_grids(pdf_plot)
        add_default_axes(pdf_plot)

        #***************************Label*************************************
        pdf_label = DataLabel(component=pdf_plot, data_point=(2.4,0.15), 
                          label_position=(15,15), padding=5,
                          label_format = 'PDF',
                          bgcolor = "transparent",
                          marker_color = "transparent",
                          marker_line_color = "transparent",
                          border_visible=False)
        pdf_plot.overlays.append(pdf_label)
        
        
#        tool = DataLabelTool(pdf_label, drag_button="right", auto_arrow_root=True)
#        pdf_label.tools.append(tool)


        container.add(pdf_plot)
        pan = PanTool(pdf_plot)
        zoom = SimpleZoom(pdf_plot, tool_mode="box", always_on=False)
        broadcaster.tools.append(pan)
        broadcaster.tools.append(zoom)
        
        
#*********************************CDF****************************

        plot = create_line_plot((self.x_array,self.pdf_array), 
                                color=tuple(COLOR_PALETTE[0]), width=2.0)
        
        plot.bgcolor = "white"
        plot.border_visible = True                

        add_default_grids(plot)
        add_default_axes(plot)


        container.add(plot)


#       # Create a pan tool and give it a reference to the plot it should
#       # manipulate, but don't attach it to the plot.  Instead, attach it to
#       # the broadcaster.
        pan = PanTool(plot)
        zoom = SimpleZoom(plot, tool_mode="box", always_on=False)
        broadcaster.tools.append(pan)
        broadcaster.tools.append(zoom)
        

    #""" PDF Plot
        # Add an axis on the right-hand side that corresponds to the second plot.
        # Note that it uses plot.value_mapper instead of plot0.value_mapper.
        pdf_plot = create_line_plot((self.x_array,self.cdf_array), 
                                color=tuple(COLOR_PALETTE[1]), width=2.0)
        pdf_plot.bgcolor = "white"
        pdf_plot.border_visible = True
     

        # Label
        cdf_text = TextBoxOverlay(text = 'CDF', alternate_position = (200,390) )
        pdf_plot.overlays.append(cdf_text)

        tool = DataLabelTool(cdf_text, drag_button="right", auto_arrow_root=True)
        cdf_text.tools.append(tool)
        
        container.add(pdf_plot)

#        vertical_axis = LabelAxis(pdf_plot, orientation='top',
#                                  title='Categories')
#        pdf_plot.underlays.append(vertical_axis)

        pdf_pan = PanTool(pdf_plot)
        pdf_zoom = SimpleZoom(pdf_plot, tool_mode="box", always_on=False)
        broadcaster.tools.append(pdf_pan)
        broadcaster.tools.append(pdf_zoom)


        axis = PlotAxis(pdf_plot, orientation="right")
        pdf_plot.underlays.append(axis)



#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

        # Add the broadcast tool to the container, instead of to an
        # individual plot
        container.tools.append(broadcaster)


        container.underlays.append(PlotLabel("CDF",
                                  component=container,
                                  font = "swiss 16",
                                  overlay_position="right"))


        legend = Legend(component=container, padding=10, align="ul")
        legend.tools.append(LegendTool(legend, drag_button="right"))
        container.overlays.append(legend)

        # Set the list of plots on the legend
        plots["pdf"] = plot
        plots["cdf"] = pdf_plot
        legend.plots = plots


#*******************************************************************************


        x = ArrayDataSource( array( [ 0.0, 0.0 ], dtype = float ) ) 
        y = ArrayDataSource( array( [ 0.0, self.mean ], dtype = float ) )
        
        
        mean_plot = create_line_plot((x,y),
                                  color=tuple(COLOR_PALETTE[2]), width=2.0)


#        vertical_plot = LinePlot(index = ArrayDataSource( array( [ 0.0, 1.0 ], dtype = float ) ),
#                        value = ArrayDataSource( array( [ 0.0, 1.0 ], dtype = float ) ),
#                        color = "green"),
#                        index_mapper = LinearMapper(range=index_mapper),
#                        value_mapper = LinearMapper(range=value_mapper))
 


        container.add(mean_plot)


       # Create a pan tool and give it a reference to the plot it should
       # manipulate, but don't attach it to the plot.  Instead, attach it to
       # the broadcaster.
        mean_pan = PanTool(mean_plot)
        mean_zoom = SimpleZoom(mean_plot, tool_mode="box", always_on=False)
        broadcaster.tools.append(mean_pan)
        broadcaster.tools.append(mean_zoom)



#**************************************************************************

        x = ArrayDataSource( array( [ 0.0, 0.0 ], dtype = float ) ) 
        y = ArrayDataSource( array( [ 0.0, 2.0 ], dtype = float ) )
        
        print "self.standard_deviation", self.standard_deviation
        
        st_plot = create_line_plot((x,y),
                                  color=tuple(COLOR_PALETTE[4]), width=2.0)

        container.add(st_plot)


       # Create a pan tool and give it a reference to the plot it should
       # manipulate, but don't attach it to the plot.  Instead, attach it to
       # the broadcaster.
        st_pan = PanTool(st_plot)
        st_zoom = SimpleZoom(st_plot, tool_mode="box", always_on=False)
        broadcaster.tools.append(st_pan)
        broadcaster.tools.append(st_zoom)



#*****************************************************************************

        # Add the title at the top
        container.overlays.append(PlotLabel("Distribution plots",
                                  component=container,
                                  font = "swiss 16",
                                  overlay_position="top"))

        # Add the traits inspector tool to the container
        container.tools.append(TraitsTool(container))
        return container

    plot = Instance(Plot)
    def _plot_default(self):
        p = Plot()
        p.tools.append(PanTool(p))
        p.overlays.append(ZoomTool(p))
        self._refresh_plot(p)
        return p

    @on_trait_change( 'loc,scale,shape,distr_type,quantile,n_segments' )
    def update_plot(self):
        p = self.plot
        p.delplot('pdf')
        p.delplot('cdf')
        self._refresh_plot(p)
    
    def _refresh_plot(self,p):
        p.datasources[ 'x' ] = ArrayDataSource( self.x_array,
                                                        sort_order = 'none' )
        
        p.datasources[ 'pdf' ] = ArrayDataSource( self.pdf_array,
                                                        sort_order = 'none' )
        p.datasources[ 'cdf' ] = ArrayDataSource( self.cdf_array,
                                                        sort_order = 'none' )
        p.plot( ('x','pdf'), name = 'pdf', color = 'blue' )
        p.plot( ('x','cdf'), name = 'cdf', color = 'red' )

    traits_view = View( HSplit( Group( Group( 
                                       Item('distr_type'), Item('scale' ), Item('shape'), Item('loc'),
                                       label =       'distribution parameters'
                                       ),
                                       Group(
                                       Item('quantile'), Item('n_segments'),
                                       label = 'plot parameters',
                                       ),
                                       Group(
                                       Item('mean',style='readonly'),
                                       Item('standard_deviation',style='readonly'),
                                       Item('skew',style='readonly'),
                                       Item('kurtosis',style='readonly'),
                                       label = 'moments',
                                       ),
                                       layout = 'tabbed',
                                       dock = 'horizontal',
                                       ),
                                       Item('plot_container',                                             
                                             editor=ComponentEditor(
                                                                    height = 400,
                                                                    width = 400,
                                                                    ), 
                                             show_label = False,
                                             resizable = True,
                                           ),
                                       Item('plot',                                             
                                             editor=ComponentEditor(
                                                                    height = 300,
                                                                    width = 300,
                                                                    ), 
                                             show_label = False,
#                                             resizable = True,
                                           ),
                                        ),
                            menubar=MenuBar(Menu(Action(name="Print data",
                                                        action="print_data"),
                                                Action(name = 'Save figure',
                                                   action='save_figure'),
                                             name="View")),
                            dock = 'tab',
                            id = 'pdistrib.Distributions',                                       
                        buttons= [OKButton,CancelButton],
                        resizable=True,
#                        width=900, height=800
                        )


if __name__ == '__main__':
    normal = PDistrib()
    normal.configure_traits()
    
    