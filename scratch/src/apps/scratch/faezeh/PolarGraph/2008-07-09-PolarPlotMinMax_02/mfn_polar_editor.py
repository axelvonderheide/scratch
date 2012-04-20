"""
Traits UI editor for WX, based on the Chaco PlotEditor in
traits.ui.wx.plot_editor.
"""


# Enthought library imports
from enthought.enable.traits.ui.wx.rgba_color_editor import \
    RGBAColorEditor
from enthought.enable.api import black_color_trait, LineStyle, ColorTrait, white_color_trait
from enthought.enable.wx_backend.api import Window
from enthought.kiva.traits.kiva_font_trait import KivaFont
from enthought.traits.api import Enum, false, Str, Range, Tuple, Float, cached_property, \
                                 Bool, Trait, Int, Any, Property, List, HasPrivateTraits, Instance, Array
from enthought.traits.ui.api import Item
from enthought.traits.ui.wx.editor import Editor
from enthought.traits.ui.wx.editor_factory import EditorFactory


# Local relative imports
from enthought.chaco.axis import PlotAxis
from enthought.chaco.plot_containers import OverlayPlotContainer

from enthought.chaco.plot_factory import create_line_plot, create_scatter_plot, \
                         add_default_grids, add_default_axes

from enthought.chaco.plot_label import PlotLabel
from enthought.chaco.scatter_markers import marker_trait

# Somewhat unorthodox...
from enthought.chaco.tools.api import PanTool, SimpleZoom, DataLabelTool
from enthought.chaco.api import create_line_plot, add_default_axes, \
                                 add_default_grids, OverlayPlotContainer, \
                                 PlotLabel, VPlotContainer, DataLabel
                                 
from enthought.chaco.chaco_plot_editor import ChacoPlotEditor
from enthought.chaco.polar_line_renderer import PolarLineRenderer
from enthought.chaco.abstract_plot_renderer import AbstractPlotRenderer
from enthought.chaco.polar_mapper import PolarMapper
from enthought.chaco.array_data_source import ArrayDataSource
from enthought.chaco.data_range_1d import DataRange1D

from numpy import frompyfunc, add, ndarray, arange, array, compress, concatenate, cos, pi, sin, transpose, zeros


#-------------------------------------------------------------------------------
#  Constants:
#-------------------------------------------------------------------------------

WindowColor = "lightgray"

#-------------------------------------------------------------------------------
#  Trait definitions:
#-------------------------------------------------------------------------------

# Range of values for an axis.
AxisRange =  Tuple( ( 0.0, 1.0, 0.01 ),
                    labels = [ 'Low', 'High', 'Step' ],
                    cols   = 3 )

# Range of axis bounds.
AxisBounds = Tuple( ( 0.0, 1.0 ),
                    labels = [ 'Min', 'Max' ],
                    cols   = 2 )

# Range for the height and width for the plot widget.
PlotSize = Range( 50, 1000, 180 )

# Range of plot line weights.
LineWeight = Range( 1, 9, 3 )

# The color editor to use for various color traits.
color_editor = RGBAColorEditor()


USE_DATA_UPDATE = 1



def create_p_plot(data, orientation='h', color='black', width=1.0,
                      dash="solid", grid="dot", value_mapper_class=PolarMapper, **kwargs):
    if (type(data) != ndarray) and (len(data) == 2):
        data = transpose(array(data))
    
    r_data, t_data = transpose(data)
    index_data= r_data*cos(t_data)
    value_data= r_data*sin(t_data)
    
    index = ArrayDataSource(index_data, sort_order='ascending')
    # Typically the value data is unsorted
    value = ArrayDataSource(value_data)

    index_range = DataRange1D()
    index_range.add(index)
    index_mapper = PolarMapper(range=index_range)
    
    value_range = DataRange1D()
    value_range.add(value)
    value_mapper = value_mapper_class(range=value_range)
    
    plot = PolarLineRenderer(index=index, value=value,
                    index_mapper = index_mapper,
                    value_mapper = value_mapper,
                    orientation = orientation,
                    color = color,
                    line_width = width,
                    line_style = dash,
                    grid_style = grid, 
                    #grid_visible= False, 
                    origin_axis_visible=True)
    return plot          


#def get_radius_min_max(radius_values):
#    r_min = radius_values.min()
#    r_max = radius_values.max()
#    return r_min, r_max
#
#def coord_trans_plt(r_value, r_min, r_max):
#    r_a = 0.3   
#    r_value_plt = r_a + (1. - r_a) * (r_value - r_min) / (r_max - r_min)
#    return r_value_plt
#
#vcoord_trans_plt = frompyfunc( coord_trans_plt, 3, 1 )
#
#
#def get_radius_values_plt(radius_values):
#    r_min, r_max = get_radius_min_max(radius_values)
#    print('r_min, r_max', r_min, r_max) 
#    radius_values_plt = vcoord_trans_plt(radius_values, r_min, r_max)
#    return radius_values_plt

#def unitcircle_fn( theta ):
#    unitcircle_plt = 0.3 
#    return unitcircle_plt
#        
#vunitcircle_fn = frompyfunc( unitcircle_fn, 1, 1 )



def get_radius_min_max(radius_arr):
    r_min = radius_arr.min()
    r_max = radius_arr.max()
    return r_min, r_max


def coord_trans_plt(r_value, r_min, r_max):
    r_a = 0.3   
    if (r_max == r_min):
        r_value_plt = r_a 
    else:
        r_value_plt = r_a + (1. - r_a) * (r_value - r_min) / (r_max - r_min)
    return r_value_plt

vcoord_trans_plt = frompyfunc( coord_trans_plt, 3, 1 )


#        radius_arr_plt = Property( Array, depends_on = 'theta_arr')
#        @cached_property
def _get_radius_arr_plt(radius_arr):
    r_min, r_max = get_radius_min_max(radius_arr)
    print 'r_min, r_max', r_min, r_max 
    radius_arr_plt = array(vcoord_trans_plt(radius_arr, r_min, r_max), dtype='float_' )
    return radius_arr_plt



def unitcircle_fn( theta_value ):
    r_a = 0.3
    unitcircle_plt = r_a 
    return unitcircle_plt
        
vunitcircle_fn = frompyfunc( unitcircle_fn, 1, 1 )

#        #unitcircle
#        unitcircle_arr = Property( Array, depends_on = 'theta_arr')
#        @cached_property
def _get_unitcircle_arr(theta_arr):
    return array( vunitcircle_fn( theta_arr ), dtype='float_' )


class MFnPolarPlotItem(Item):
    """ A Traits UI Item for a Chaco plot, for use in Traits UI Views.
    
    NOTE: ComponentEditor is preferred over this class, as it is more
    flexible.
    """
    # Name of the trait that references the index data source.
    index = Str
    # Name of the trait that references the value data source.
    value_list = List( Str )
    # Title of the plot (overlaid on the plot container).
    title = Str("Polar Plot Editor")

    # Bounds of the x-axis, used if **x_auto** is False.
    x_bounds = AxisBounds
    # Set the x-axis bounds automatically?
    x_auto = Bool(True)
    # Bounds of the y-axis, used if **y_auto** is False.
    y_bounds = AxisBounds
    # Set the y-axis bounds automatically?
    y_auto = Bool(True)

    # The orientation of the index axis.
    orientation = Enum("h", "v")

    # If these are None, then the index/value trait names are used

    # Label of the x-axis; if None, the **index** name is used.
    x_label = Trait(None, None, Str)
    
    # Name of the trait on the object containing the label of the x-axis.
    # This takes precedence over **x_label**.
    x_label_trait = Trait(None, None, Str)
    
    # Font for the label of the x-axis.
    x_label_font = KivaFont("modern 10")
    # Color of the label of the x-axis.
    x_label_color = black_color_trait
    # Label of the y-axis; if None, the **value** name is used.
    y_label = Trait(None, None, Str)
    # Name of the trait on the object containing the label of the y-axis.
    # This takes precedence over **y_label**.
    y_label_trait = Trait(None, None, Str)
    # Font for the label of the y-axis.
    y_label_font = KivaFont("modern 10")
    # Color of the label of the y-axis.
    y_label_color = black_color_trait

    # General plot properties

    # Foreground olor of the plot.
    color = ColorTrait("green")
    # Background color of the plot.
    bgcolor = white_color_trait
    # Background color of the plot (deprecated).
    bg_color = Property   # backwards compatibility; deprecated
    # Color of the background padding.
    padding_bg_color = ColorTrait(WindowColor)

    # Border properties

    # Width of the plot border
    border_width = Int(1)
    # Is the border visible?
    border_visible = false
    # Line style of the border.
    border_dash = LineStyle
    # Color of the border.
    border_color = black_color_trait

    # The type of the plot.
    type = Enum("line", "scatter")
    # The type of the plot as a string.
    type_trait = Str

    # plot-specific properties.  These might not apply to all plot types.

    # Type of marker (for plots that use markers).
    marker = marker_trait
    # Size of marker (for plots that use markers).
    marker_size = Int(4)
    # Marker outline color (for plots that user markers).
    outline_color = black_color_trait


    def __init__(self, index, value_list, type="line", **traits):
        self.index = index
        self.value_list = value_list
        self.type = type
        self.name = "Plot"
        super(Item, self).__init__(**traits)

        self.editor = PolarEditorFactory()

        self.editor.plotitem = self

        return


    def _set_bg_color(self, val):
        self.bgcolor = val

    def _get_bg_color(self):
        return self.bgcolor


class PolarEditorFactory ( EditorFactory ):
    """ Editor factory for plot editors.
    """
    #---------------------------------------------------------------------------
    #  Trait definitions:
    #---------------------------------------------------------------------------

    # Width of the plot editor.
    width    = PlotSize
    # Height of the plot editor.
    height   = PlotSize
    # The ChacoPlotItem associated with this factory.
    plotitem = Any


    #---------------------------------------------------------------------------
    #  'Editor' factory methods:
    #---------------------------------------------------------------------------

    def simple_editor ( self, ui, object, name, description, parent ):
        return MFnPolarPlotEditor( parent,
                                 factory     = self,
                                 ui          = ui,
                                 object      = object,
                                 name        = name,
                                 description = description )

    def text_editor ( self, ui, object, name, description, parent ):
        return MFnPolarPlotEditor( parent,
                                 factory     = self,
                                 ui          = ui,
                                 object      = object,
                                 name        = name,
                                 description = description )

    def readonly_editor ( self, ui, object, name, description, parent ):
        return MFnPolarPlotEditor( parent,
                                 factory     = self,
                                 ui          = ui,
                                 object      = object,
                                 name        = name,
                                 description = description )


class MFnPolarPlotEditor ( Editor ):
    """ Traits UI editor for displaying trait values in a Chaco plot.
    """
    # MFnPolarPlotEditorToolbar associated with the editor:
    toolbar = Any

 

    def init ( self, parent ):
        """ Finishes initializing the editor by creating the underlying toolkit
            widget.
        """
        factory = self.factory
        plotitem = factory.plotitem

        container = OverlayPlotContainer(padding = 50, fill_padding = True,
                                         bgcolor = plotitem.padding_bg_color,
                                         use_backbuffer=True)

        if plotitem.title != '':
            container.overlays.append(PlotLabel(plotitem.title, component=container,
                                                overlay_position="top"))

        self._container = container
        window = Window(parent, component = container)
        self.control = control = window.control
        control.SetSize((factory.width, factory.height))
        object = self.object

        # Attach listeners to the object's traits appropriately so we can
        # update the plot when they change.  For the _update_axis_grids()
        # callback, we have to wrap it in a lambda to keep traits from
        # inferring the calling convention based on introspecting the argument
        # list.
        for name in [plotitem.index] + plotitem.value_list:
            print "setting on trait change for ", name
            object.on_trait_change( self._update_data, name)
        object.on_trait_change(self.update_editor, plotitem.type_trait)
        return

    #---------------------------------------------------------------------------
    #  Disposes of the contents of an editor:
    #---------------------------------------------------------------------------

    def dispose(self):
        """ Disposes of the contents of the editor.
        """
        object = self.object
        plotitem = self.factory.plotitem

        if USE_DATA_UPDATE == 1:
            for name in (plotitem.index, plotitem.value_list):
                object.on_trait_change( self._update_data, name, remove = True )
        for name in (plotitem.type_trait,):
            object.on_trait_change( self._update_editor, name, remove = True )
        self._destroy_plot()
        super(MFnPolarPlotEditor, self).dispose()

    def _destroy_plot(self):
        if self._container and self._plot:
            plot = self._plot
            del plot.index._data
            del plot.index._cached_mask
            del plot.value._data
            del plot.value._cached_mask
            self._container.remove(plot)
            self._plot = None
            plot.index = None
            plot.value = None

#  Destroy the second Plot if two curves (functions) are to be displayed
#            plot = self._plot2
#            del plot.index._data
#            del plot.index._cached_mask
#            del plot.value._data
#            del plot.value._cached_mask
#            self._container.remove(plot)
#            self._plot = None
#            plot.index = None
#            plot.value = None

        return

 
    #---------------------------------------------------------------------------
    #  Finishes initializing the editor by creating the underlying toolkit
    #  widget:
    #---------------------------------------------------------------------------

    def update_editor(self):
        """ Updates the editor when the object trait changes externally to the
            editor.
        """

        factory  = self.factory
        plotitem = factory.plotitem

        # Remove the old plot
        if self._plot is not None:
            self._destroy_plot()

        try:
            theta_arr = getattr(self.object, plotitem.index)
            radius_arr = getattr(self.object, plotitem.value_list[0])
#            plotrange_min = getattr(self.object, plotitem.value_list[1])
#            plotrange_max = getattr(self.object, plotitem.value_list[2])
            
        except:
            self._container.request_redraw()
            return

        if plotitem.type_trait != "":
            plot_type = getattr(self.object, plotitem.type_trait)
        else:
            plot_type = plotitem.type

        if plotitem.x_auto == True:
            index_bounds = None
        else:
            index_bounds = plotitem.x_bounds

        if plotitem.y_auto == True:
            value_bounds = None
        else:
            value_bounds = plotitem.y_bounds



        unitcircle_arr = _get_unitcircle_arr(theta_arr)

        print 'unitcircle_arr', unitcircle_arr

        plot2 = self._create_polar_plot( plotitem, (unitcircle_arr, theta_arr),
                                         orientation = plotitem.orientation, color='red', 
                                         width=5.0 )
        
        
        self._plot2 = plot2
        self._container.add(plot2)
        self._container.request_redraw()

       
        radius_arr_plt = _get_radius_arr_plt(radius_arr)

        print 'radius_arr_plt', radius_arr_plt

        plot = self._create_polar_plot( plotitem, (radius_arr_plt, theta_arr),
                                        orientation = plotitem.orientation,
                                        width=5.0, index_sort="ascending" )

        self._set_basic_properties(plot, plotitem)       

        self._plot = plot
        self._container.add(plot)
        self._container.request_redraw()


        
        
        # Add some tools
#        plot2.tools.append(PanTool(plot2))
#        zoom = SimpleZoom(plot2, tool_mode="box", always_on=False)
#        plot2.overlays.append(zoom)
#        
#        
#        label = DataLabel(component=plot2, data_point=(unitcircle_values, theta_values),
#                          label_position="top left", padding=40,
#                          bgcolor = "lightgray",
#                          border_visible=False)
#        plot2.overlays.append(label)
#        tool = DataLabelTool(label, drag_button="right", auto_arrow_root=True)
#        label.tools.append(tool)

        

 


    def _update_data(self):
        """ Updates the editor when the object trait changes externally to the
            editor.
        """
        self.update_editor()

    #---------------------------------------------------------------------------
    #  Finishes initializing the editor by creating the underlying toolkit
    #  widget:
    #---------------------------------------------------------------------------
        """ Finishes initializing the editor by creating the underlying toolkit
            widget.
        """
        
    #---------------------------------------------------------------------------
    #  Creates the table editing tool bar:
    #---------------------------------------------------------------------------

    def _create_toolbar ( self, parent, sizer ):
        """ Creates the table editing toolbar.
        """
        factory = self.factory
        if not factory.show_toolbar:
            return

    def _create_polar_plot(self, plotitem, values, **kwargs):
        plot = create_p_plot(values, **kwargs)
        return plot          


    def _create_bar_plot(self, plotitem, values, **kwargs):
        plot = create_bar_plot(values, **kwargs)
        return plot            


    def _set_basic_properties(self, plot, plotitem):
        for attr in ("color", "bgcolor", "border_visible", "border_width",
                     "border_dash", "border_color"):
            setattr(plot, attr, getattr(plotitem, attr))
        return

    def _create_line_plot(self, plotitem, values, **kwargs):
        plot = create_line_plot(values, **kwargs)
        return plot

    def _create_scatter_plot(self, plotitem, values, **kwargs):
        plot = create_scatter_plot(values, **kwargs)
        for attr in ("marker", "marker_size", "outline_color"):
            setattr(plot, attr, getattr(plotitem, attr))
        return plot

    def _add_axis_grids(self, new_plot, plotitem):
        value_axis, index_axis = add_default_axes(new_plot, 
                                    orientation=plotitem.orientation)
        add_default_grids(new_plot)
        new_plot.tools.append(PanTool(new_plot))
        zoom = SimpleZoom(component=new_plot, tool_mode="box", always_on=False)
        new_plot.overlays.append(zoom)

        # Update the titles and labels
        self._update_axis_grids(new_plot, plotitem)

    def _update_axis_grids(self, plot=None, plotitem=None):
        if plot is None:
            plot = self._plot
        if plotitem is None:
            plotitem = self.factory.plotitem

        if plotitem.x_label_trait is not None:
            htitle = getattr(self.object, plotitem.x_label_trait)
        elif plotitem.x_label is not None:
            htitle = plotitem.x_label
        else:
            htitle = plotitem.index

        if plotitem.y_label_trait is not None:
            htitle = getattr(self.object, plotitem.y_label_trait)
        elif plotitem.y_label is not None:
            vtitle = plotitem.y_label
        else:
            vtitle = plotitem.value_list

        if plotitem.orientation == "v":
            htitle, vtitle = vtitle, htitle
        plot.x_axis.title = htitle
        plot.y_axis.title = vtitle
        
        # This is sort of crappy.. since we are using BaseXYPlots and not
        # Plot/DataViews, we actually can't easily get references to the plot's
        # index and value axes.  So we have to search through the underlays for
        # PlotAxis instances whose ranges match the index and value ranges.
        for axis in plot.underlays + plot.overlays:
            if isinstance(axis, PlotAxis) and axis.mapper.range is plot.index_range:
                axis.title_font = plotitem.x_label_font
                axis.title_color = plotitem.x_label_color

        for axis in plot.underlays + plot.overlays:
            if isinstance(axis, PlotAxis) and axis.mapper.range is plot.value_range:
                axis.title_font = plotitem.y_label_font
                axis.title_color = plotitem.y_label_color

        
        plot.request_redraw()
        return


class MFnPolarPlotEditorToolbar ( HasPrivateTraits ):
    """ Toolbar displayed in table editors.
    """
    
     # The table editor that this is the toolbar for:
    editor = Instance( MFnPolarPlotEditor )

    # The toolbar control:
    control = Any

    #---------------------------------------------------------------------------
    #  Initializes the toolbar for a specified window:
    #---------------------------------------------------------------------------

    def __init__ ( self, parent = None, **traits ):
        super( MFnPolarPlotEditorToolbar, self ).__init__( **traits )
        factory = self.editor.factory
        actions = []
        
        if factory.sortable and (not factory.sort_model):
            actions.append( self.no_sort )
            
        if self.editor.in_row_mode:
            if factory.reorderable:
                actions.append( self.move_up )
                actions.append( self.move_down )
                
            if factory.search is not None:
                actions.append( self.search )
            
        if factory.editable:
            if (factory.row_factory is not None) and (not factory.auto_add):
                actions.append( self.add )
            if ((factory.deletable != False) and 
                (factory.selection_mode in ( 'row', 'rows' ))):
                actions.append( self.delete )
                
        if factory.configurable:
            actions.append( self.prefs )
            
        if len( actions ) > 0:
            self.control.SetBackgroundColour( parent.GetBackgroundColour() )

            # fixme: Why do we have to explictly set the size of the toolbar?
            #        Is there some method that needs to be called to do the
            #        layout?
            self.control.SetSize( wx.Size( 23 * len( actions ), 16 ) )

    #---------------------------------------------------------------------------
    #  PyFace/Traits menu/toolbar controller interface:
    #---------------------------------------------------------------------------

    #---------------------------------------------------------------------------
    #  Adds a menu item to the menu bar being constructed:
    #---------------------------------------------------------------------------

    def add_to_menu ( self, menu_item ):
        """ Adds a menu item to the menu bar being constructed.
        """
        pass

    #---------------------------------------------------------------------------
    #  Adds a tool bar item to the tool bar being constructed:
    #---------------------------------------------------------------------------

    def add_to_toolbar ( self, toolbar_item ):
        """ Adds a toolbar item to the too bar being constructed.
        """
        pass

    #---------------------------------------------------------------------------
    #  Returns whether the menu action should be defined in the user interface:
    #---------------------------------------------------------------------------

    def can_add_to_menu ( self, action ):
        """ Returns whether the action should be defined in the user interface.
        """
        return True

    #---------------------------------------------------------------------------
    #  Returns whether the toolbar action should be defined in the user
    #  interface:
    #---------------------------------------------------------------------------

    def can_add_to_toolbar ( self, action ):
        """ Returns whether the toolbar action should be defined in the user
            interface.
        """
        return True

    #---------------------------------------------------------------------------
    #  Performs the action described by a specified Action object:
    #---------------------------------------------------------------------------

    def perform ( self, action, action_event = None ):
        """ Performs the action described by a specified Action object.
        """
        getattr( self.editor, action.action )()



class PolarLineRenderer(AbstractPlotRenderer):
    """ A renderer for polar line plots.
    """
    #------------------------------------------------------------------------
    # Appearance-related traits
    #------------------------------------------------------------------------

    # The color of the origin axis.
    origin_axis_color_ = (0,0,0,1)
    # The width of the origin axis.
    origin_axis_width = 1.0
    # The origin axis is visible.
    origin_axis_visible=True
    # The grid is visible.
    grid_visible= True
    # The orientation of the plot is horizontal; for any other value, it is 
    # transposed 
    orientation = 'h'
    # The color of the line.
    color = black_color_trait
    # The width of the line.
    line_width = Float(1.0)
    # The style of the line.
    line_style = LineStyle("solid")
    # The style of the grid lines.
    grid_style= LineStyle("dot")

    def _gather_points(self):
        """
        Collects the data points that are within the plot bounds and caches them
        """
        # This is just a stub for now.  We should really find the lines only
        # inside the screen range here.

        x = self.index.get_data()
        y = self.value.get_data()
        rad= min(self.width/2.0,self.height/2.0)
        sx = x*rad+ self.x + self.width/2.0
        sy = y*rad+ self.y + self.height/2.0

        points = transpose(array((sx,sy)))
        self._cached_data_pts = points
        self._cache_valid = True
        return

    def _data_changed(self):
        self._cache_valid = False
        return

    def _update_mappers(self):
        #Dunno if there is anything else to do here
        self._cache_valid = False

    def _render(self, gc, points):
        """ Actually draw the plot.
        """
        gc.save_state()

        gc.set_antialias(True)
        self._draw_default_axes(gc)
        self._draw_default_grid(gc)
        if len(points)>0:
            gc.clip_to_rect(self.x, self.y, self.width, self.height)
            gc.set_stroke_color(self.color_)
            gc.set_line_width(self.line_width)
            gc.set_line_dash(self.line_style_)

            gc.begin_path()
            gc.lines(points)
            gc.stroke_path()

        gc.restore_state()

        return

    def map_screen(self, data_array):
        """ Maps an array of data points into screen space and returns it as
        an array. 
        
        Implements the AbstractPlotRenderer interface.
        """

        if len(data_array) == 0:
            return []
        elif len(data_array) == 1:
            xtmp, ytmp = transpose(data_array)
            x_ary = xtmp
            y_ary = ytmp
        else:
            x_ary, y_ary = transpose(data_array)

        sx = self.index_mapper.map_screen(x_ary)
        sy = self.value_mapper.map_screen(y_ary)

        if self.orientation == 'h':
            return transpose(array((sx, sy)))
        else:
            return transpose(array((sy, sx)))

    def map_data(self, screen_pt):
        """ Maps a screen space point into the "index" space of the plot.
        
        Implements the AbstractPlotRenderer interface.
        """
        if self.orientation == 'h':
            x, y = screen_pt
        else:
            y,x = screen_pt
        return array((self.index_mapper.map_data(x),
                      self.value_mapper.map_data(y)))


    def _downsample(self):
        return self.map_screen(self._cached_data_pts)

    def _draw_plot(self, *args, **kw):
        """ Draws the 'plot' layer.
        """
        # Simple compatibility with new-style rendering loop
        return self._draw_component(*args, **kw)


    def _draw_component(self, gc, view_bounds=None, mode='normal'):
        """ Renders the component. 
        """
        self._gather_points()
        self._render(gc, self._cached_data_pts)

    def _bounds_changed(self, old, new):
        super(PolarLineRenderer, self)._bounds_changed(old, new)
        self._update_mappers()

    def _bounds_items_changed(self, event):
        super(PolarLineRenderer, self)._bounds_items_changed(event)
        self._update_mappers()

    def _draw_default_axes(self, gc):
        if not self.origin_axis_visible:
            return
        gc.save_state()
        gc.set_stroke_color(self.origin_axis_color_)
        gc.set_line_width(self.origin_axis_width)
        gc.set_line_dash(self.grid_style_)
        x_data,y_data= transpose(self._cached_data_pts)
        x_center=self.x + self.width/2.0
        y_center=self.y + self.height/2.0

        # number of divisions used to divide the cirlce radially
        # (equal to number of axes to be plotted)
        n_axes = 2
        for theta in range(n_axes*2):
                r= min(self.width/2.0,self.height/2.0)
                x= r*cos(theta*pi/n_axes) + x_center
                y= r*sin(theta*pi/n_axes) + y_center
                data_pts= array([[x_center,y_center],[x,y]])
                start,end = data_pts
                gc.move_to(int(start[0]), int(start[1]))
                gc.line_to(int(end[0]), int(end[1]))
                gc.stroke_path()
  
#        for theta in range(12):
#                r= min(self.width/2.0,self.height/2.0)
#                x= r*cos(theta*pi/6) + x_center
#                y= r*sin(theta*pi/6) + y_center
#                data_pts= array([[x_center,y_center],[x,y]])
#                start,end = data_pts
#                gc.move_to(int(start[0]), int(start[1]))
#                gc.line_to(int(end[0]), int(end[1]))
#                gc.stroke_path()
  
        gc.restore_state()
        return

    def _draw_default_grid(self,gc):
        if not self.grid_visible:
            return
        gc.save_state()
        gc.set_stroke_color(self.origin_axis_color_)
        gc.set_line_width(self.origin_axis_width)
        gc.set_line_dash(self.grid_style_)
        x_data,y_data= transpose(self._cached_data_pts)
        x_center=self.x + self.width/2.0
        y_center=self.y + self.height/2.0
#        for r_part in range(5):
#            rad= min(self.width/2.0,self.height/2.0)
#            r= rad*(0.3 + 0.7 * r_part/4) 
#            gc.move_to(self.x,self.y)
#            gc.arc(x_center,y_center,r,0,2*pi)
#            gc.stroke_path()

        
       
       # Only two grids  
#        gc.move_to(self.x,self.y)
#
#        rad_one = min(self.width/2.0,self.height/2.0)
#        rad_unitcircle_plt = 0.3
#        rad_0_75 = (rad_unitcircle_plt + (1.0 - rad_unitcircle_plt) * 0.75 ) * rad_one
#        rad_0_50 = (rad_unitcircle_plt + (1.0 - rad_unitcircle_plt) * 0.50 ) * rad_one
#        rad_0_25 = (rad_unitcircle_plt + (1.0 - rad_unitcircle_plt) * 0.25 ) * rad_one
#        rad_zero = (rad_unitcircle_plt + (1.0 - rad_unitcircle_plt) * 0.00 ) * rad_one

        # number of divisions used to divide the cirlce tangentially
        # (equal to number of gridlines to be plotted)
#        ndiv = 4
#        for i in range(ndiv):
#                rad_i= (0.3 + 0.7 * (1/ndiv)*(i+1) ) * rad_one
#                gc.arc(x_center,y_center,rad_i,0,2*pi)
#                gc.stroke_path()

        rad_one = min(self.width/2.0,self.height/2.0)
        rad_unitcircle_plt = 0.3   
        n_grids = 5
        for i in range(n_grids):
            rad_i         = (0.3 + 0.7*  1/n_grids  * (i+1) ) * rad_one     
            gc.move_to(self.x,self.y)
            gc.arc(x_center,y_center,rad_i,0,2*pi)
            gc.stroke_path()

 
#        #internal grid (value = 0)
#        gc.arc(x_center,y_center,rad_zero,0,2*pi)
#        gc.stroke_path()
#
#        #internal grid (value = 0.25)
#        gc.arc(x_center,y_center,rad_0_25,0,2*pi)
#        gc.stroke_path()
#
#       #internal grid (value = 0.5)
#        gc.arc(x_center,y_center,rad_0_50,0,2*pi)
#        gc.stroke_path()
#
#       #internal grid (value = 0.75)
#        gc.arc(x_center,y_center,rad_0_75,0,2*pi)
#        gc.stroke_path()
#
#
#        #external grid (Value = 1)
#        gc.arc(x_center,y_center,rad_one,0,2*pi)
#        gc.stroke_path()

        gc.restore_state()
        return



# EOF
