"""
Traits UI editor for WX, based on the Chaco PlotEditor in
traits.ui.wx.plot_editor.
"""

# Enthought library imports
from enthought.enable.api import black_color_trait, LineStyle, ColorTrait, white_color_trait
from enthought.enable.wx_backend.api import Window

from enthought.traits.api import false, Str, Range, Float, Bool, Int, Any, List
from enthought.traits.ui.api import Item
from enthought.traits.ui.wx.editor import Editor
from enthought.traits.ui.wx.editor_factory import EditorFactory

# Local relative imports
from enthought.chaco.plot_containers import OverlayPlotContainer
from enthought.chaco.plot_label import PlotLabel

# Somewhat unorthodox...
from enthought.chaco.tools.api import SimpleZoom, DataLabelTool
from enthought.chaco.api import OverlayPlotContainer, PlotLabel, Label, DataLabel
from enthought.chaco.polar_line_renderer import PolarLineRenderer
from enthought.chaco.abstract_plot_renderer import AbstractPlotRenderer
from enthought.chaco.polar_mapper import PolarMapper
from enthought.chaco.array_data_source import ArrayDataSource
from enthought.chaco.data_range_1d import DataRange1D

from numpy import frompyfunc, add, ndarray, arange, array, compress, concatenate, cos, pi, sin, transpose, zeros


#Constant variable
WindowColor = "lightgray"

# Range for the height and width for the plot widget.
PlotSize = Range( 50, 1000, 180 )


def coord_trans_plt(r_value, frac_noplot, plotrange_min, plotrange_max):
    r_a = frac_noplot   
    if r_value < plotrange_min :
        r_value_plt = r_a
    elif r_value > plotrange_max :
        r_value_plt = 1.0
    else:
        r_value_plt = r_a + (1. - r_a) * (r_value - plotrange_min) / (plotrange_max - plotrange_min)
    return r_value_plt




def _get_radius_arr_plt(radius_arr, frac_noplot, plotrange_min, plotrange_max):
    vcoord_trans_plt = frompyfunc( coord_trans_plt, 4, 1 )
    radius_arr_plt = array(vcoord_trans_plt(radius_arr, frac_noplot, plotrange_min, plotrange_max), dtype='float_' )
    return radius_arr_plt


def unitcircle_fn( theta_value, frac_noplot, plotrange_min, plotrange_max):
    r_a = coord_trans_plt(0.0, frac_noplot, plotrange_min, plotrange_max)
    unitcircle_plt = r_a 
    return unitcircle_plt
        
def _get_unitcircle_arr(theta_arr, frac_noplot, plotrange_min, plotrange_max):
    vunitcircle_fn = frompyfunc( unitcircle_fn, 4, 1 )    
    return array( vunitcircle_fn( theta_arr, frac_noplot, plotrange_min, plotrange_max), dtype='float_' )


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


    # Foreground olor of the plot.
    color = ColorTrait("green")
    # Background color of the plot.
    bgcolor = white_color_trait
    # Color of the background padding.
    padding_bg_color = ColorTrait(WindowColor)

    orientation = "h"


    # Width of the plot border
    border_width = Int(1)
    # Is the border visible?
    border_visible = false
    # Color of the border.
    border_color = black_color_trait

    # The type of the plot as a string.
    type_trait = Str


    def __init__(self, index, value_list, type="line", **traits):
        self.index = index
        self.value_list = value_list
#        self.type = type
        self.name = "Plot"

        super(Item, self).__init__(**traits)
        self.editor = MFnPolarEditorFactory()
        self.editor.plotitem = self

        return


    def _set_bg_color(self, val):
        self.bgcolor = val

    def _get_bg_color(self):
        return self.bgcolor


class MFnPolarEditorFactory ( EditorFactory ):

    # Width of the plot editor.
    width    = PlotSize
    # Height of the plot editor.
    height   = PlotSize
    # The ChacoPlotItem associated with this factory.
    plotitem = Any

    def simple_editor ( self, ui, object, name, description, parent ):
        return MFnPolarPlotEditor( parent,
                                 factory     = self,
                                 ui          = ui,
                                 object      = object,
                                 name        = name,
                                 description = description )
        

class MFnPolarPlotEditor ( Editor ):
    """ Traits UI editor for displaying trait values in a Chaco plot.
    """

    # flag that indicates whether the unitcircle is to be plotted
    # (i.e. if zero is in the plotrange)
    plot_unitcircle_flag = Bool(False)
    
    def validate(self, plotrange_min, plotrange_max):
        if plotrange_min >= plotrange_max:
            print "### Error ###: invalid plot ranges. plotrange_max must be greater than plotrange_min."
            return False
        else:
            return True    

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
            theta_arr     = getattr(self.object, plotitem.index)
            radius_arr    = getattr(self.object, plotitem.value_list[0])
            plotrange_min = getattr(self.object, plotitem.value_list[1])
            plotrange_max = getattr(self.object, plotitem.value_list[2])
            frac_noplot   = getattr(self.object, plotitem.value_list[3])
            
        except:
            self._container.request_redraw()
            return

        # check plausibility
        if not ( (array([plotrange_min]) <= radius_arr).all() and (radius_arr <= array([plotrange_max])).all() ):
            print "Note: some value out of range! Compare plotranges with indicated values for Radius_min and Radius_max "



        if self.validate(plotrange_min, plotrange_max): 
            # plot unitcircle if zero is in plotrange
            unitcircle_arr = _get_unitcircle_arr(theta_arr, frac_noplot, plotrange_min, plotrange_max)
    
            bool_arr_frac_noplot = (unitcircle_arr == array([frac_noplot]))
            bool_arr_1           = (unitcircle_arr == array([1.]))
            
            if ( plotrange_min == 0. or plotrange_max == 0.) or  \
               ( bool_arr_frac_noplot.all() == False and bool_arr_1.all() == False ):
                self.plot_unitcircle_flag = True
                plot_unitcircle = self._create_polar_plot( plotitem, (unitcircle_arr, theta_arr),
                                             orientation = plotitem.orientation, 
                                             color='red', 
                                             frac_noplot = frac_noplot,
                                             width=5.0 )
                
                self._plot_unitcircle = plot_unitcircle
                self._container.add(plot_unitcircle)
                self._container.request_redraw()
                
   
#            radius_axeslabel_min = array([ frac_noplot ])
#            theta_axeslabel_min  = array([      0      ])
#            
#            plot_axeslabel_min = self._create_polar_plot( plotitem, (radius_axeslabel_min, theta_axeslabel_min),
#                                         orientation = plotitem.orientation, 
#                                         color='red', 
#                                         frac_noplot = frac_noplot,
#                                         width=5.0 )
#            
            ### label
#            label = Label(component=plot_axeslabel_min, 
#                              label_position="top left", padding=40,
#                              bgcolor = "lightgray",
#                              border_visible=False)
#            
#            pltlblright = PlotLabel("1", component=plot_axeslabel_min, vjustify="center", overlay_position="right")
#            pltlblcenter = PlotLabel("0", component=plot_axeslabel_min, vjustify="center") 
#            
#            self._theta_axeslabel_min = theta_axeslabel_min
#            
#            self._container.overlays.append(pltlblright)
#            self._container.overlays.append(pltlblcenter)
#            self._container.add(plot_axeslabel_min)
#            self._container.request_redraw()
            
   
    
            # function plot
            radius_arr_plt = _get_radius_arr_plt(radius_arr, frac_noplot, plotrange_min, plotrange_max)
    
            plot_radiusfn = self._create_polar_plot( plotitem, (radius_arr_plt, theta_arr),
                                            orientation = plotitem.orientation,
                                            frac_noplot = frac_noplot,
                                            width=5.0, 
                                            index_sort="ascending" )
    
            self._set_basic_properties(plot_radiusfn, plotitem)       
    
            self._plot_radiusfn = plot_radiusfn
            self._container.add(plot_radiusfn)
            self._container.request_redraw()
        else: 
            self._plot_radiusfn = None
    
    
            
            
            
#        plot_unitcircle.tools.append(PanTool(plot_unitcircle))
#        zoom = SimpleZoom(plot_unitcircle, tool_mode="box", always_on=False)
#        plot_unitcircle.overlays.append(zoom)
#            
            
#        label = DataLabel(component=plot_unitcircle, data_point=(unitcircle_values, theta_values),
#                              label_position="top left", padding=40,
#                              bgcolor = "lightgray",
#                              border_visible=False)
#        plot_unitcircle.overlays.append(label)
        
#        tool = DataLabelTool(label, drag_button="right", auto_arrow_root=True)
#        label.tools.append(tool)
        

 


    def _update_data(self):
        """ Updates the editor when the object trait changes externally to the
            editor.
        """
        self.update_editor()

    #---------------------------------------------------------------------------
    #  Creates the table editing tool bar:
    #---------------------------------------------------------------------------

    def _create_polar_plot(self, plotitem, data, 
                           orientation='h', 
                           color='black', 
                           width=1.0,
                           dash="solid", 
                           grid="dot", 
                           value_mapper_class = PolarMapper,
                           frac_noplot = 0.3,
                           **kwargs):
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
        
        plot = MFnPolarLineRenderer(index=index, value=value,
                                 index_mapper = index_mapper,
                                 value_mapper = value_mapper,
                                 orientation = orientation,
                                 color = color,
                                 line_width = width,
                                 line_style = dash,
                                 grid_style = grid,
                                 frac_noplot = frac_noplot, 
                                 origin_axis_visible=True)
        return plot          

    def _set_basic_properties(self, plot, plotitem):
        for attr in ("color", "bgcolor", "border_visible", "border_width",
                     "border_color"):
                #"border_dash", "border_color"):
            setattr(plot, attr, getattr(plotitem, attr))
        return

#class MFnPolarPlotEditorToolbar ( HasPrivateTraits ):
#    """ Toolbar displayed in table editors.
#    """
#    
#     # The table editor that this is the toolbar for:
#    editor = Instance( MFnPolarPlotEditor )
#
#    # The toolbar control:
#    control = Any
#
#    #---------------------------------------------------------------------------
#    #  Initializes the toolbar for a specified window:
#    #---------------------------------------------------------------------------
#
#    def __init__ ( self, parent = None, **traits ):
#        super( MFnPolarPlotEditorToolbar, self ).__init__( **traits )
#        factory = self.editor.factory
#        actions = []
#        
#        if factory.sortable and (not factory.sort_model):
#            actions.append( self.no_sort )
#            
#        if self.editor.in_row_mode:
#            if factory.reorderable:
#                actions.append( self.move_up )
#                actions.append( self.move_down )
#                
#            if factory.search is not None:
#                actions.append( self.search )
#            
#        if factory.editable:
#            if (factory.row_factory is not None) and (not factory.auto_add):
#                actions.append( self.add )
#            if ((factory.deletable != False) and 
#                (factory.selection_mode in ( 'row', 'rows' ))):
#                actions.append( self.delete )
#                
#        if factory.configurable:
#            actions.append( self.prefs )
#            
#        if len( actions ) > 0:
#            self.control.SetBackgroundColour( parent.GetBackgroundColour() )
#
#            # fixme: Why do we have to explictly set the size of the toolbar?
#            #        Is there some method that needs to be called to do the
#            #        layout?
#            self.control.SetSize( wx.Size( 23 * len( actions ), 16 ) )
#
#    #---------------------------------------------------------------------------
#    #  PyFace/Traits menu/toolbar controller interface:
#    #---------------------------------------------------------------------------
#
#    #---------------------------------------------------------------------------
#    #  Adds a menu item to the menu bar being constructed:
#    #---------------------------------------------------------------------------
#
#    def add_to_menu ( self, menu_item ):
#        """ Adds a menu item to the menu bar being constructed.
#        """
#        pass
#
#    #---------------------------------------------------------------------------
#    #  Adds a tool bar item to the tool bar being constructed:
#    #---------------------------------------------------------------------------
#
#    def add_to_toolbar ( self, toolbar_item ):
#        """ Adds a toolbar item to the too bar being constructed.
#        """
#        pass
#
#    #---------------------------------------------------------------------------
#    #  Returns whether the menu action should be defined in the user interface:
#    #---------------------------------------------------------------------------
#
#    def can_add_to_menu ( self, action ):
#        """ Returns whether the action should be defined in the user interface.
#        """
#        return True
#
#    #---------------------------------------------------------------------------
#    #  Returns whether the toolbar action should be defined in the user
#    #  interface:
#    #---------------------------------------------------------------------------
#
#    def can_add_to_toolbar ( self, action ):
#        """ Returns whether the toolbar action should be defined in the user
#            interface.
#        """
#        return True
#
#    #---------------------------------------------------------------------------
#    #  Performs the action described by a specified Action object:
#    #---------------------------------------------------------------------------
#
#    def perform ( self, action, action_event = None ):
#        """ Performs the action described by a specified Action object.
#        """
#        getattr( self.editor, action.action )()



class MFnPolarLineRenderer(AbstractPlotRenderer):
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
    
    frac_noplot = Float(0.3)

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
        # (equals the number of axes to be plotted)
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

        rad_one = min(self.width/2.0,self.height/2.0)
        rad_unitcircle_plt = 0.3   

        # number of divisions used to divide the cirlce tangentially
        ndivs_grid = 4
        for i in range(ndivs_grid+1):
            plotrange_max = 1.0
            plotrange_min = 0.0
            rad_i = plotrange_min + (plotrange_max-plotrange_min) / ndivs_grid * i
            rad_i_plt =  coord_trans_plt(rad_i, self.frac_noplot, plotrange_min, plotrange_max) * rad_one
            gc.move_to(self.x,self.y)
            gc.arc(x_center, y_center, rad_i_plt, 0, 2*pi)
            gc.stroke_path()


        gc.restore_state()
        return



# EOF
