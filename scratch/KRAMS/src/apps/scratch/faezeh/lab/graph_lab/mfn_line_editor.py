# Enthought library imports
from enthought.enable.traits.ui.wx.rgba_color_editor import \
    RGBAColorEditor
from enthought.enable.api import black_color_trait, LineStyle, ColorTrait, white_color_trait
from enthought.enable.wx_backend.api import Window
from enthought.kiva.traits.kiva_font_trait import KivaFont
from enthought.traits.api import Enum, false, Str, Range, Tuple, \
                                 Bool, Trait, Int, Any, Property, Instance, HasPrivateTraits
from enthought.traits.ui.api import Item, UI
from enthought.traits.ui.menu import Action, ToolBar, Menu
from enthought.traits.ui.wx.editor import Editor
from enthought.traits.ui.wx.editor_factory import EditorFactory
from enthought.traits.ui.wx.helper import traits_ui_panel


from enthought.chaco.api import create_line_plot, add_default_axes, \
                                 add_default_grids, OverlayPlotContainer, \
                                 PlotLabel, VPlotContainer, \
                                 create_scatter_plot, Legend, PlotComponent, PlotGraphicsContext
from enthought.chaco.tools.api import PanTool, RectZoomTool, SimpleZoom, \
                                       LegendTool, TraitsTool, DragZoom

# Somewhat unorthodox...
from enthought.chaco.tools.api import PanTool, SimpleZoom
# Chaco imports
from enthought.chaco.chaco_plot_editor import ChacoPlotEditor, \
                                                ChacoPlotItem

#from enthought.pyface.dock.core \
#    import DockWindow, DockSizer, DockSection, DockRegion, DockControl

#from enthought.traits.ui.dock_window_theme \
#    import DockWindowTheme

from numpy import savetxt, hstack, vstack

import wx

from enthought.pyface.api import FileDialog, OK

from enthought.pyface.image_resource \
    import ImageResource

import csv

USE_DATA_UPDATE = 1
WILDCARD = "Saved plots (*.eps)|*.eps|"\
           "All files (*.*)|*.*"

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

    # The DockWindow graphical theme:
    #dock_theme = Instance( DockWindowTheme )
 
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

    # TableEditorToolbar associated with the editor:
    toolbar = Any

    # The Traits UI associated with the table editor toolbar:
    toolbar_ui = Instance( UI )


    # The Traits UI associated with the function editor toolbar:
    #toolbar_ui = Instance( UI )

    #---------------------------------------------------------------------------
    #  Finishes initializing the editor by creating the underlying toolkit
    #  widget:
    #---------------------------------------------------------------------------

    def update_editor( self):
        super(MFnWTPlotEditor,self).update_editor()
        self._plot.y_axis.tick_label_formatter = lambda val: "%g" % val
        self._plot.x_axis.tick_label_formatter = lambda val: "%g" % val

    #---------------------------------------------------------------------------
    #  Finishes initializing the editor by creating the underlying toolkit
    #  widget:
    #---------------------------------------------------------------------------
        """ Finishes initializing the editor by creating the underlying toolkit
            widget.
        """
        
    def init ( self, parent ):

        factory = self.factory
        plotitem = factory.plotitem

        panel = traits_ui_panel( parent, -1 )

        super(MFnWTPlotEditor,self).init(panel)

        plot_control = self.control
        self.control = panel

        sizer        = wx.BoxSizer( wx.VERTICAL )
        self._create_toolbar( panel, sizer )
        
        sizer.Add( plot_control, 1, wx.ALIGN_RIGHT | wx.EXPAND )

        panel.SetSizer( sizer )
 
        return

        
    #---------------------------------------------------------------------------
    #  Creates the table editing tool bar:
    #---------------------------------------------------------------------------

    def _create_toolbar ( self, parent, sizer ):
        """ Creates the table editing toolbar.
        """
        factory = self.factory
            
        toolbar = MFnWTPlotEditorToolbar( parent = parent, editor = self )
        tb_sizer = wx.BoxSizer( wx.HORIZONTAL )

        self.toolbar = toolbar
        tb_sizer.Add( toolbar.control, 0 )
        tb_sizer.Add( ( 1, 1 ), 1, wx.EXPAND )

        sizer.Add( wx.StaticLine( parent, -1, style = wx.LI_HORIZONTAL ), 0,
                   wx.EXPAND | wx.BOTTOM, 5 )
        sizer.Add( tb_sizer, 0, wx.ALIGN_RIGHT | wx.EXPAND )
        sizer.Add( wx.StaticLine( parent, -1, style = wx.LI_HORIZONTAL ), 0,
                   wx.EXPAND | wx.BOTTOM, 5 )

    #---------------------------------------------------------------------------
    #  Handles the user requesting that columns not be sorted:
    #---------------------------------------------------------------------------

    def on_savedata ( self ):
        """ Handles the user requesting that the data of the function is to be saved.
        """
        import os
        dlg = FileDialog(parent = self.control, 
                         title = 'Export function data',
                         default_directory=os.getcwd(),
                         default_filename="", wildcard='*.csv',
                         action='save as')
        if dlg.open() == OK:
            path = dlg.path

            print "Saving data to", path, "..."
            try:

                factory  = self.factory
                plotitem = factory.plotitem
                x_values = getattr(self.object, plotitem.index)
                y_values = getattr(self.object, plotitem.value)
                savetxt( path, vstack( (x_values,y_values) ).transpose() )
            except:
                print "Error saving!"
                raise
            print "Plot saved."
        return

    def on_savefig ( self ):
        """ Handles the user requesting that the image of the function is to be saved.
        """
        import os
        dlg = FileDialog(parent = self.control, 
                         title = 'Save as image',
                         default_directory=os.getcwd(),
                         default_filename="", wildcard=WILDCARD,
                         action='save as')
        if dlg.open() == OK:
            path = dlg.path

            print "Saving plot to", path, "..."
            try:

                # Now we create a canvas of the appropriate size and ask it to render
                # our component.  (If we wanted to display this plot in a window, we
                # would not need to create the graphics context ourselves; it would be
                # created for us by the window.)
                self._plot.bounds = [500,300]
                self._plot.padding = 50
                plot_gc = PlotGraphicsContext(self._plot.outer_bounds)
                print self._plot.outer_bounds
                plot_gc.render_component(self._plot)

                # Finally, we tell the graphics context to save itself to disk as an image.
                plot_gc.save(path)
                
            except:
                print "Error saving!"
                raise
            print "Plot saved."
        return

#-------------------------------------------------------------------------------
#  'MFnWTPlotEditorToolbar' class:
#-------------------------------------------------------------------------------

class MFnWTPlotEditorToolbar ( HasPrivateTraits ):
    """ Toolbar displayed in table editors.
    """
    
    #---------------------------------------------------------------------------
    #  Trait definitions:
    #---------------------------------------------------------------------------

    # Do not sort columns:
    save_data = Instance( Action,
                        { 'name':    'Save as data',
                         'tooltip': 'Save the function values',
                         'action':  'on_savedata',
                         'enabled': True,
                         'image':   ImageResource( 'table_no_sort.png' ) } )

     # Move current object up one row:
    save_fig = Instance( Action,
                        { 'name':    'Save as fig',
                         'tooltip': 'Save as figure',
                         'action':  'on_savefig',
                         'enabled': True,
                         'image':   ImageResource( 'table_move_down.png' ) } )

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

        actions = [ self.save_data, self.save_fig ]
        toolbar = ToolBar( image_size      = ( 16, 16 ),
                           show_tool_names = False,
                           show_divider    = False,
                           *actions )
        self.control = toolbar.create_tool_bar( parent, self )
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

