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
#LineWeight = Range( 1, 9, 3 )

# The color editor to use for various color traits.
#color_editor = RGBAColorEditor()


USE_DATA_UPDATE = 1

WILDCARD = "Saved plots (*.eps)|*.eps|"\
           "All files (*.*)|*.*"



class MFnMultiLinePlotItem(Item):
    """ A Traits UI Item for a Chaco plot, for use in Traits UI Views.
    
    NOTE: ComponentEditor is preferred over this class, as it is more
    flexible.
    """
    # Name of the trait that references the index data source.
    index = Str
    # Name of the trait that references the value data source.
    value_list = List( Str )
    # Title of the plot (overlaid on the plot container).
    title = Str("Multiline Plot Editor")

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

        self.editor = MultiLineEditorFactory()

        self.editor.plotitem = self

        return


    def _set_bg_color(self, val):
        self.bgcolor = val

    def _get_bg_color(self):
        return self.bgcolor


class MultiLineEditorFactory ( EditorFactory ):
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
        return MFnMultiLinePlotEditor( parent,
                                 factory     = self,
                                 ui          = ui,
                                 object      = object,
                                 name        = name,
                                 description = description )

    def text_editor ( self, ui, object, name, description, parent ):
        return MFnMultiLinePlotEditor( parent,
                                 factory     = self,
                                 ui          = ui,
                                 object      = object,
                                 name        = name,
                                 description = description )

    def readonly_editor ( self, ui, object, name, description, parent ):
        return MFnMultiLinePlotEditor( parent,
                                 factory     = self,
                                 ui          = ui,
                                 object      = object,
                                 name        = name,
                                 description = description )


class MFnMultiLinePlotEditor ( Editor ):
    """ Traits UI editor for displaying trait values in a Chaco plot.
    """
    # MFnMultiLinePlotEditorToolbar associated with the editor:
    toolbar = Any

    frac_noplot = Float(0.3)
    
    plot_unitcircle_flag = Bool(False)
    
#    if (self.unitcircle_arr.all() != frac_noplot) and (self.unitcircle_arr.all() != 1.0):
#        plot_unitcircle_flag = True
        

    def init ( self, parent ):
        """ Finishes initializing the editor by creating the underlying toolkit
            widget.
        """
        factory = self.factory
        plotitem = factory.plotitem

        container = OverlayPlotContainer(padding = 50, fill_padding = True,
                                         bgcolor = plotitem.padding_bg_color,
                                         use_backbuffer=True)

#        if plotitem.title != '':
#            container.overlays.append(PlotLabel(plotitem.title, component=container,
#                                                overlay_position="top"))

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
        super(MFnMultiLinePlotEditor, self).dispose()

    def _destroy_plot(self):
        plot_fn_list = [self._plot_radiusfn]
        if self.plot_unitcircle_flag:
            plot_fn_list += [self._plot_unitcircle]
            self.plot_unitcircle_flag = False
        for plot_fn in plot_fn_list:
            if self._container and plot_fn:
                plot = plot_fn
                del plot.index._data
                del plot.index._cached_mask
                del plot.value._data
                del plot.value._cached_mask
                self._container.remove(plot)
                plot_fn = None
                plot.index = None
                plot.value = None
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
        if self._plot_radiusfn is not None:
            self._destroy_plot()

        try:
            x_data     = getattr(self.object, plotitem.index)
            y_data     = getattr(self.object, plotitem.value_list[0]) 
            
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


        plot_y_data = self._create_line_plot( plotitem, (y_data, x_data),
                                        orientation = plotitem.orientation,
                              #          frac_noplot = frac_noplot,
                                        width=5.0, 
                                        index_sort="ascending" )

        self._set_basic_properties(plot_y_data, plotitem)       

        self._plot_y_data = plot_y_data
        self._container.add(plot_y_data)
        self._container.request_redraw()


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
            plot = self._plot_radiusfn
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


class MFnMultiLinePlotEditorToolbar ( HasPrivateTraits ):
    """ Toolbar displayed in table editors.
    """
    
     # The table editor that this is the toolbar for:
    editor = Instance( MFnMultiLinePlotEditor )

    # The toolbar control:
    control = Any

    #---------------------------------------------------------------------------
    #  Initializes the toolbar for a specified window:
    #---------------------------------------------------------------------------

    def __init__ ( self, parent = None, **traits ):
        super( MFnMultiLinePlotEditorToolbar, self ).__init__( **traits )
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




# EOF
