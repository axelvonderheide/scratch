# Enthought library imports
#from enthought.enable.traits.ui.wx.rgba_color_editor import \
#    RGBAColorEditor
from enthought.enable.api import black_color_trait, LineStyle, ColorTrait, white_color_trait
from enthought.enable.wx_backend.api import Window
from enthought.kiva.traits.kiva_font_trait import KivaFont
from enthought.traits.api import Enum, false, Str, Range, Tuple, \
                                 Bool, Trait, Int, Any, Property, Instance, HasPrivateTraits
from enthought.traits.ui.api import Item
from enthought.traits.ui.menu import Action
from enthought.traits.ui.wx.editor import Editor
from enthought.traits.ui.wx.editor_factory import EditorFactory


from enthought.chaco.api import create_line_plot, add_default_axes, \
                                 create_polar_plot, \
                                 add_default_grids, OverlayPlotContainer, \
                                 PlotLabel, VPlotContainer, \
                                 create_scatter_plot, Legend, PlotComponent
from enthought.chaco.tools.api import PanTool, RectZoomTool, SimpleZoom, \
                                       LegendTool, TraitsTool, DragZoom

# Somewhat unorthodox...
from enthought.chaco.tools.api import PanTool, SimpleZoom
# Chaco imports
from enthought.chaco.chaco_plot_editor import ChacoPlotEditor, \
                                                ChacoPlotItem

import wx

#from helper \
#    import open_fbi, Orientation, traits_ui_panel

# Range for the height and width for the plot widget.
PlotSize = Range( 50, 1000, 180 )

class MFnWTPlotItem( ChacoPlotItem ):
    def __init__(self, index, value, type="line", **traits):
        super(MFnWTPlotItem, self).__init__(index,value,type,**traits)
        self.editor = MFnWTEditorFactory()
        self.editor.plotitem = self
        return

class MFnWTEditorFactory ( EditorFactory ):
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
        return MFnWTPlotEditor( parent,
                                 factory     = self,
                                 ui          = ui,
                                 object      = object,
                                 name        = name,
                                 description = description )

    def text_editor ( self, ui, object, name, description, parent ):
        return MFnWTPlotEditor( parent,
                                 factory     = self,
                                 ui          = ui,
                                 object      = object,
                                 name        = name,
                                 description = description )

    def readonly_editor ( self, ui, object, name, description, parent ):
        return MFnWTPlotEditor( parent,
                                 factory     = self,
                                 ui          = ui,
                                 object      = object,
                                 name        = name,
                                 description = description )

class MFnWTPlotEditor ( ChacoPlotEditor ):
    """ Traits UI editor for displaying trait values in a MFnWT.
    """

    # MFnWTPlotEditorToolbar associated with the editor:
    toolbar = Any

    # The Traits UI associated with the function editor toolbar:
    #toolbar_ui = Instance( UI )

    #---------------------------------------------------------------------------
    #  Finishes initializing the editor by creating the underlying toolkit
    #  widget:
    #---------------------------------------------------------------------------

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
        for name in (plotitem.index, plotitem.value):
            object.on_trait_change( self._update_data, name)
        object.on_trait_change(self.update_editor, plotitem.type_trait)
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
            x_values = getattr(self.object, plotitem.index)
            y_values = getattr(self.object, plotitem.value)
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

        # Class-level attribute mapping different plot_type strings to methods for
        # creating different types of plots
        plot_creator_map = {"polar": self._create_polar_plot,
                            "line": self._create_line_plot,
                            "scatter": self._create_scatter_plot }

        if plot_type == 'polar':
            print
            print x_values
            print
            print y_values
            plot = self._create_polar_plot( plotitem, (y_values, x_values),
                                            orientation = plotitem.orientation,
                                            width=5 )
        elif plot_type in plot_creator_map.keys():
            plot = plot_creator_map[plot_type](plotitem, (x_values, y_values),
                                                index_bounds = index_bounds,
                                                value_bounds = value_bounds,
                                                orientation = plotitem.orientation)
            self._add_axis_grids(plot, plotitem)
            self._plot.y_axis.tick_label_formatter = lambda val: "%g" % val
            self._plot.x_axis.tick_label_formatter = lambda val: "%g" % val

        else:
            raise RuntimeError, "Unknown plot type '%s' in ChacoPlotEditor." % plot_type

        self._set_basic_properties(plot, plotitem)

        self._plot = plot
        self._container.add(plot)
        self._container.request_redraw()

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
        plot = create_polar_plot(values, **kwargs)
        return plot            
#-------------------------------------------------------------------------------
#  'MFnWTPlotEditorToolbar' class:
#-------------------------------------------------------------------------------

class MFnWTPlotEditorToolbar ( HasPrivateTraits ):
    """ Toolbar displayed in table editors.
    """
    
    #---------------------------------------------------------------------------
    #  Trait definitions:
    #---------------------------------------------------------------------------

##     # Do not sort columns:
##     no_sort = Instance( Action,
##                         { 'name':    'No Sorting',
##                           'tooltip': 'Do not sort columns',
##                           'action':  'on_no_sort',
##                           'enabled': False,
##                           'image':   ImageResource( 'table_no_sort.png' ) } )

##     # Move current object up one row:
##     move_up = Instance( Action,
##                         { 'name':    'Move Up',
##                           'tooltip': 'Move current item up one row',
##                           'action':  'on_move_up',
##                           'enabled': False,
##                           'image':   ImageResource( 'table_move_up.png' ) } )

##     # Move current object down one row:
##     move_down = Instance( Action,
##                           { 'name':    'Move Down',
##                             'tooltip': 'Move current item down one row',
##                             'action':  'on_move_down',
##                             'enabled': False,
##                             'image':   ImageResource( 'table_move_down.png' ) })

##     # Search the table:
##     search = Instance( Action,
##                        { 'name':    'Search',
##                          'tooltip': 'Search table',
##                          'action':  'on_search',
##                          'image':   ImageResource( 'table_search.png' ) } )

##     # Add a row:
##     add = Instance( Action,
##                     { 'name':    'Add',
##                       'tooltip': 'Insert new item',
##                       'action':  'on_add',
##                       'image':   ImageResource( 'table_add.png' ) } )

##     # Delete selected row:
##     delete = Instance( Action,
##                        { 'name':    'Delete',
##                          'tooltip': 'Delete current item',
##                          'action':  'on_delete',
##                          'image':   ImageResource( 'table_delete.png' ) } )

##     # Edit the user preferences:
##     prefs = Instance( Action,
##                       { 'name':    'Preferences',
##                         'tooltip': 'Set user preferences for table',
##                         'action':  'on_prefs',
##                         'image':   ImageResource( 'table_prefs.png' ) } )

    # The table editor that this is the toolbar for:
    editor = Instance( MFnWTPlotEditor )

    # The toolbar control:
    control = Any

    #---------------------------------------------------------------------------
    #  Initializes the toolbar for a specified window:
    #---------------------------------------------------------------------------

    def __init__ ( self, parent = None, **traits ):
        super( MFnWTPlotEditorToolbar, self ).__init__( **traits )
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

