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

    # a button to reset the graph
    refresh_view = Button('Refresh View')
    
    #The reset function for the above button  
    def _refresh_view_fired(self):
        self.update_container()
        
    #the reset function for menu bar. It resets the graph
    def reset_view(self):
        self.update_container()

    #Exit the form
    def exit_file(self):
        sys.exit(0)

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
        self._refresh_container( container )
        return container

    @on_trait_change( 'loc,scale,shape,distr_type,quantile,n_segments' )
    def update_container(self):
        c = self.plot_container
        c.remove( *c.components )
        self._refresh_container(c)

    def reset_view(self):
        self.update_container()
    
    def _refresh_container(self,container):
        """ Plot some distribution functions """
        plots = {}
        broadcaster = BroadcasterTool()
    
        index = ArrayDataSource(self.x_array)
        value = ArrayDataSource(self.pdf_array, sort_order="none")

        index_range = DataRange1D()
        index_range.add(index)
        index_mapper = LinearMapper(range=index_range)
    
        value_range = DataRange1D(  low_setting = 0.0  )
        value_range.add(value)
        value_mapper = LinearMapper(range=value_range)

        """Plot probability distribution function(pdf) with marking the area under the function  """ 
        plot_pdf = FilledLinePlot(index = index, value = value,
                                    index_mapper = index_mapper,
                                    value_mapper = value_mapper,
                                    edge_color = tuple(COLOR_PALETTE[0]),
                                    face_color = "paleturquoise",
                                    border_visible = True)
        
        """define the grid, axes and title of the vertical grid """
        add_default_grids(plot_pdf)
        add_default_axes(plot_pdf, vtitle="PDF") 

        """create a label for the pdf and append it to the plot_pdf """
        label_pdf = DataLabel(component=plot_pdf, data_point=(2.4,0.13), 
                          label_position=(15,15), padding=5,
                          label_format = 'PDF',
                          bgcolor = "transparent",
                          marker_color = "transparent",
                          marker_line_color = "transparent",
                          arrow_color = tuple(COLOR_PALETTE[0]),
                          border_visible=False)
        plot_pdf.overlays.append(label_pdf)
        container.add(plot_pdf)
        
        """create a label for the x coordinate """
        container.overlays.append(PlotLabel("X",
                                  component=container,
                                  font = "swiss 16",
                                  overlay_position="bottom"))

        pan = PanTool(plot_pdf)
#        zoom = SimpleZoom(plot_pdf, tool_mode="box", always_on=False)

        broadcaster.tools.append(pan)
 #       broadcaster.tools.append(zoom)
        #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
        #""" Mean Plot
        x = ArrayDataSource( array( [ self.mean, self.mean ], dtype = float ) ) 
        y = ArrayDataSource( array( [ 0.0, max( self.pdf_array ) ], dtype = float ) )

        """ Plot the mean value"""
        plot_mean = LinePlot(index = x, value = y,
                        color = "pink",
                        index_mapper = index_mapper,
                        value_mapper = value_mapper )
        container.add(plot_mean)
        
        # Create a pan tool and give it a reference to the plot it should
        # manipulate, but don't attach it to the plot.  Instead, attach it to
        # the broadcaster.
        mean_pan = PanTool(plot_mean)
#        mean_zoom = SimpleZoom(plot_mean, tool_mode="box", always_on=False)
        broadcaster.tools.append(mean_pan)
 #       broadcaster.tools.append(mean_zoom)
        #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
        # Add an axis on the right-hand side that corresponds to the second plot.
        # Note that it uses plot.value_mapper instead of plot0.value_mapper.
        """ Plot cdf and its label """
        plot_cdf = create_line_plot((self.x_array,self.cdf_array), 
                                color=tuple(COLOR_PALETTE[1]), width=2.0)
        plot_cdf.bgcolor = "white"
        plot_cdf.border_visible = True

        label_cdf = DataLabel(component=plot_cdf, data_point=(2.4,0.9), 
                          label_position=(-35,20), padding=5,
                          label_format = 'CDF',
                          bgcolor = "transparent",
                          marker_color = "transparent",
                          marker_line_color = "transparent",
                          arrow_color = tuple(COLOR_PALETTE[1]),
                          border_visible=False)
        plot_cdf.overlays.append(label_cdf)

        container.add(plot_cdf)

        pan1 = PanTool(plot_cdf)
#        zoom1 = SimpleZoom(plot_cdf, tool_mode="box", always_on=False)
        broadcaster.tools.append(pan1)
#        broadcaster.tools.append(zoom1)

        axis = PlotAxis(plot_cdf, title="CDF", orientation="right")
        plot_cdf.underlays.append(axis)
        #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
        # Add the broadcast tool to the container, instead of to an
        # individual plot
        container.tools.append(broadcaster)
        ##**************************************************************************
        """Plot standard deviation """
        stdev_low   = self.mean - self.standard_deviation
        stdev_high  = self.mean + self.standard_deviation
        stdev_low_height = self.distr.pdf( stdev_low )
        stdev_high_height = self.distr.pdf( stdev_high )
        
        x = ArrayDataSource( array( [ stdev_low, stdev_low, stdev_high, stdev_high ], dtype = float ) ) 
        y = ArrayDataSource( array( [ 0, stdev_low_height, stdev_high_height, 0. ], dtype = float ) )
              
        plot_st = PolygonPlot(index = x, value = y,
                        edge_color = "purple",
                        index_mapper = index_mapper,
                        value_mapper = value_mapper )
        container.add(plot_st)
        ##*****************************************************************************
        """create the legend for the mean """
        legend = Legend(component=container, padding=10, align="ul")
        legend.tools.append(LegendTool(legend, drag_button="right"))
        container.overlays.append(legend)

        # Set the list of plots on the legend
        plots["Mean"] = plot_mean
        legend.plots = plots

        # Add the title at the top
        container.overlays.append(PlotLabel("probability distribution plots",
                                  component=container,
                                  font = "swiss 16",
                                  overlay_position="top"))
        # Add the traits inspector tool to the container
        container.tools.append(TraitsTool(container))

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

    traits_view = View( HSplit( 
                               Group( Group( 
                                       Item('distr_type'), Item('scale' ), Item('shape'), Item('loc'),
                                       Item('refresh_view',style='simple', label='Reset View'), 
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
                                        ),
                            menubar=MenuBar(Menu(Action(name ="Print data",action="print_data"),
                                                 Action(name ="R&eset", action="reset_view"),
                                                 Action(name ="E&xit",action="exit_file"),
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