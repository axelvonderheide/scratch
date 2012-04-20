#-------------------------------------------------------------------------------
#
# Copyright (c) 2009, IMB, RWTH Aachen.
# All rights reserved.
#
# This software is provided without warranty under the terms of the BSD
# license included in simvisage/LICENSE.txt and may be redistributed only
# under the conditions described in the aforementioned license.  The license
# is also available online at http://www.simvisage.com/licenses/BSD.txt
#
# Created on May 25, 2009 by Rostislav Chudoba
#
#-------------------------------------------------------------------------------

# Enthought library imports
from enthought.traits.api import \
    Enum, false, Str, Range, Tuple, Bool, Trait, Int, Any, Property, Instance, HasPrivateTraits
from enthought.traits.ui.api import  Item, UI
from enthought.traits.ui.menu import Action, ToolBar, Menu
from enthought.traits.ui.wx.editor import Editor
from enthought.traits.ui.basic_editor_factory import BasicEditorFactory

from enthought.enable.api import black_color_trait, LineStyle, ColorTrait, white_color_trait, Window, Component
from enthought.kiva.traits.kiva_font_trait import KivaFont
from enthought.kiva.backend_image import GraphicsContext

from enthought.chaco.example_support import COLOR_PALETTE
from enthought.chaco.api import \
    Plot, ArrayPlotData, create_line_plot, add_default_axes, PlotLabel, VPlotContainer,\
    add_default_grids, OverlayPlotContainer, PlotAxis, Legend, ArrayPlotData, ArrayDataSource, \
    create_scatter_plot, Legend, PlotComponent, PlotGraphicsContext, DataRange1D, LinePlot, LinearMapper, \
    Plot, ArrayPlotData
from enthought.chaco.tools.api import \
                                PanTool, RectZoomTool, SimpleZoom, \
                                LegendTool, TraitsTool, DragZoom, \
                                BroadcasterTool, ZoomTool, LegendTool

from enthought.pyface.api import FileDialog, OK, ImageResource
from enthought.pyface.image_resource import ImageResource

from mfn_plot_adapter import MFnPlotAdapter

from scipy.special import jn
import wx, csv #, sys

from numpy import *

USE_DATA_UPDATE = 1
WILDCARD = "Saved plots (*.png)|*.png|"\
           "All files (*.*)|*.*"

# Range for the height and width for the plot widget.
PlotSize = Range( 50, 1000, 180 )

class _MFnChacoEditor ( Editor ):
    """ Traits UI editor for displaying trait values in a MFnLine.
    """
  
    # TableEditorToolbar associated with the editor:
    toolbar = Any
    # The Traits UI associated with the function editor toolbar:
    #  toolbar_ui = Instance( UI )

    # adjustable parameters
    adapter = Instance( MFnPlotAdapter )
    
    splot = Instance(Plot)
    line_plot = Instance(LinePlot)
    
    #---------------------------------------------------------------------------
    #  Finishes initializing the editor by creating the underlying toolkit
    #  widget:
    #---------------------------------------------------------------------------

    plot_container = Instance( OverlayPlotContainer )
    def _plot_container_default( self ):
        container = OverlayPlotContainer( padding = 70, fill_padding = False,
                                  bgcolor = "white", use_backbuffer=True)

        return container
    
    #---------------------------------------------------------------------------
    #  Finishes initializing the editor by creating the underlying toolkit
    #  widget:
    #---------------------------------------------------------------------------
        """ Finishes initializing the editor by creating the underlying toolkit
            widget.
        """
        
    def init ( self, parent ):

        factory = self.factory
        self.adapter = factory.adapter

        self.control = self._create_canvas(parent)

        # Register the update listener
        #
        self.value.on_trait_change( self.update_editor, 'data_changed' )

    def update_editor( self):
        c = self.plot_container
        c.remove( *c.components )
        self._refresh_container()        

    def _create_canvas( self, parent ):
        '''Create canvas for chaco plots
        '''
        panel = wx.Panel( parent, -1, style=wx.CLIP_CHILDREN )
        
        container_panel = Window( panel, component = self.plot_container )
        
        sizer = wx.BoxSizer( wx.VERTICAL )
        sizer.Add( wx.StaticLine( parent, -1, style = wx.LI_HORIZONTAL ), 0,
                   wx.EXPAND | wx.BOTTOM, 5 )
        sizer.Add( self._create_toolbar( panel ), 0, wx.EXPAND )
        sizer.Add( container_panel.control, 1, wx.EXPAND )

        sizer.Add( wx.StaticLine( parent, -1, style = wx.LI_HORIZONTAL ), 0,
                   wx.EXPAND | wx.BOTTOM, 5 )

        panel.SetSizer(sizer)

        return panel

    def _refresh_container( self ):
        ''' rebuild the container for the current data
        '''
        broadcaster = BroadcasterTool()
    
#        mfn_line = self.value
#        ydata = transpose(mfn_line.ydata)
        
        adapter = self.adapter
        if adapter.var_x != '':
            # Get the x-label text from the object's trait var_x          
            label_x = getattr( self.object, adapter.var_x )
        else:
            # Get the x-label from the adapter          
            label_x = adapter.label_x
            
        if adapter.var_y  != '':
            label_y = getattr( self.object, adapter.var_y )
        else:
            label_y = adapter.label_y

        i=0          # for colors
#        plots = {}   # for legend
#        index = ArrayDataSource(mfn_line.xdata)

#****************************************** Massivbau Table
        x = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.296, 0.30, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.371 ]
        index = ArrayDataSource(x)
        index_range = DataRange1D()
        index_range.add(index)
        index_mapper = LinearMapper(range=index_range)
            
        value_range = DataRange1D(  low_setting = 0.0  )
        
        


            
        vector1 = [0.030, 0.044, 0.055, 0.066, 0.076, 0.086, 0.097, 0.107, 0.118, 0.131, 0.145, 0.159, 0.173, 0.188, 0.202, 0.217, 0.232, 0.248, 0.250, 0.264, 0.280, 0.296, 0.312, 0.329, 0.346, 0.364, 0.382, 0.400, 0.419, 0.438, 0.450, 0.458, 0.478, 0.499, 0.520, 0.542, 0.565, 0.589, 0.614, 0.617]
        y1 = ArrayDataSource(vector1, sort_order="none")
        value_range.add(y1)
        value_mapper = LinearMapper(range=value_range)
        
        self.line_plot = LinePlot(index = index, value = y1,
                                        index_mapper = index_mapper,
                                        value_mapper = value_mapper,
                                        color = tuple(COLOR_PALETTE[i]),
                                        width = 20,
                                        edge_color = 'blue',
                                        border_visible = False)
    
        add_default_grids(self.line_plot)
        add_default_axes(self.line_plot, vtitle="x/d", htitle="x")
    
        self.plot_container.add(self.line_plot)
        

#        omega = [0.0101, 0.0203, 0.0306, 0.0410, 0.0515, 0.0621, 0.0728, 0.0836, 0.0946, 0.1057, 0.1170, 0.1285, 0.1401, 0.1518, 0.1638, 0.1759, 0.1882, 0.2007, 0.2024, 0.2134, 0.2263, 0.2395, 0.2529, 0.2665, 0.2804, 0.2946, 0.3091, 0.3239, 0.3391, 0.3546, 0.3643, 0.3706, 0.3869, 0.4038, 0.4211, 0.4391, 0.4576, 0.4768, 0.4968, 0.4994 ]
#        y_omega = ArrayDataSource(omega, sort_order="none")
#        value_range.add(y_omega)
#        value_mapper = LinearMapper(range=value_range)
#        
#        self.line_plot = LinePlot(index = index, value = y_omega,
#                                        index_mapper = index_mapper,
#                                        value_mapper = value_mapper,
#                                        color = tuple(COLOR_PALETTE[i]),
#                                        width = 20,
#                                        edge_color = 'blue',
#                                        border_visible = False)
#    
#        add_default_grids(self.line_plot)
#        add_default_axes(self.line_plot, vtitle="omega", htitle="x")
#    
#        self.plot_container.add(self.line_plot)
#        
#        
#        vector2 = [0.990, 0.985, 0.980, 0.976, 0.971, 0.967, 0.962, 0.956, 0.951, 0.946, 0.940, 0.934, 0.928, 0.922, 0.916, 0.910, 0.903, 0.897, 0.896, 0.890, 0.884, 0.877, 0.870, 0.863, 0.856, 0.849, 0.841, 0.834, 0.826, 0.818, 0.813, 0.810, 0.801, 0.793, 0.784, 0.774, 0.765, 0.755, 0.745, 0.743]
#        y2 = ArrayDataSource(vector2, sort_order="none")
#        value_range.add(y2)
#        value_mapper = LinearMapper(range=value_range)
#        
#        self.line_plot = LinePlot(index = index, value = y2,
#                                        index_mapper = index_mapper,
#                                        value_mapper = value_mapper,
#                                        color = tuple(COLOR_PALETTE[i]),
#                                        width = 20,
#                                        edge_color = 'blue',
#                                        border_visible = False)
#        add_default_grids(self.line_plot)
#        add_default_axes(self.line_plot, vtitle="z/d", htitle="x")
#        self.plot_container.add(self.line_plot)
#
#
#        vector3 = [-0.77, -1.15, -1.46, -1.76, -2.06, -2.37, -2.68, -3.01, -3.35, -3.5, -3.5, -3.5, -3.5, -3.5, -3.5, -3.5, -3.5, -3.5, -3.5, -3.5, -3.5, -3.5, -3.5, -3.5, -3.5, -3.5, -3.5, -3.5, -3.5, -3.5, -3.5, -3.5, -3.5, -3.5, -3.5, -3.5, -3.5, -3.5, -3.5, -3.5]
#        y3 = ArrayDataSource(vector3, sort_order="none")
#        value_range.add(y3)
#        value_mapper = LinearMapper(range=value_range)
#        
#        self.line_plot = LinePlot(index = index, value = y3,
#                                        index_mapper = index_mapper,
#                                        value_mapper = value_mapper,
#                                        color = tuple(COLOR_PALETTE[i]),
#                                        width = 20,
#                                        edge_color = 'blue',
#                                        border_visible = False)
#    
#        add_default_grids(self.line_plot)
#        add_default_axes(self.line_plot, vtitle="epsilon_c2", htitle="x")
#    
#        self.plot_container.add(self.line_plot)
#
#
#        vector4 = [25, 25, 25, 25, 25, 25, 25, 25, 25, 23.29, 20.71, 18.55, 16.73, 15.16, 13.80, 12.61, 11.55, 10.62, 10.50, 9.78, 9.02, 8.33, 7.71, 7.13, 6.60, 6.12, 5.67, 5.25, 4.86, 4.49, 4.28, 4.15, 3.82, 3.52, 3.23, 2.95, 2.69, 2.44, 2.20, 2.174]
#        y4 = ArrayDataSource(vector4, sort_order="none")
#        value_range.add(y4)
#        value_mapper = LinearMapper(range=value_range)
#        
#        self.line_plot = LinePlot(index = index, value = y4,
#                                        index_mapper = index_mapper,
#                                        value_mapper = value_mapper,
#                                        color = tuple(COLOR_PALETTE[i]),
#                                        width = 20,
#                                        edge_color = 'blue',
#                                        border_visible = False)
#    
#        add_default_grids(self.line_plot)
#        add_default_axes(self.line_plot, vtitle="epsilon_s1", htitle="x")
#    
#        self.plot_container.add(self.line_plot)
#
#
#        vector5 = 434.8 * ones(25)#[0.030, 0.044, 0.055, 0.066, 0.076, 0.876, 0.097, 0.107, 0.118, 0.131, 0.145, 0.159, 0.173, 0.188, 0.202, 0.217, 0.232, 0.248, 0.250, 0.264, 0.280, 0.296, 0.312, 0.329, 0.346, 0.364, 0.382, 0.400, 0.419, 0.438, 0.450, 0.458, 0.478, 0.499, 0.520, 0.542, 0.565, 0.589, 0.614, 0.617]
#        y5 = ArrayDataSource(vector5, sort_order="none")
#        value_range.add(y5)
#        value_mapper = LinearMapper(range=value_range)
#        
#        self.line_plot = LinePlot(index = index, value = y5,
#                                        index_mapper = index_mapper,
#                                        value_mapper = value_mapper,
#                                        color = tuple(COLOR_PALETTE[i]),
#                                        width = 20,
#                                        edge_color = 'blue',
#                                        border_visible = False)
#        add_default_grids(self.line_plot)
#        add_default_axes(self.line_plot, vtitle="gamma_sd", htitle="x")
#        self.plot_container.add(self.line_plot)
#
#
#        vector6 = [0.337, 0.464, 0.553, 0.622, 0.676, 0.718, 0.751, 0.778, 0.801, 0.810, 0.810, 0.810, 0.810, 0.810, 0.810, 0.810, 0.810, 0.810, 0.810, 0.810, 0.810, 0.810, 0.810, 0.810, 0.810, 0.810, 0.810, 0.810, 0.810, 0.810, 0.810, 0.810, 0.810, 0.810, 0.810, 0.810, 0.810, 0.810, 0.810, 0.810]
#        y6 = ArrayDataSource(vector6, sort_order="none")
#        value_range.add(y6)
#        value_mapper = LinearMapper(range=value_range)
#        
#        self.line_plot = LinePlot(index = index, value = y6,
#                                        index_mapper = index_mapper,
#                                        value_mapper = value_mapper,
#                                        color = tuple(COLOR_PALETTE[i]),
#                                        width = 20,
#                                        edge_color = 'blue',
#                                        border_visible = False)
#    
#        add_default_grids(self.line_plot)
#        add_default_axes(self.line_plot, vtitle="alpha_R", htitle="x")
#    
#        self.plot_container.add(self.line_plot)
#
#        vector7 = [0.346, 0.353, 0.360, 0.368, 0.377, 0.387, 0.396, 0.405, 0.413, 0.416, 0.416, 0.416, 0.416, 0.416, 0.416, 0.416, 0.416, 0.416, 0.416, 0.416, 0.416, 0.416, 0.416, 0.416, 0.416, 0.416, 0.416, 0.416, 0.416, 0.416, 0.416, 0.416, 0.416, 0.416, 0.416, 0.416, 0.416, 0.416, 0.416, 0.416]
#        y7 = ArrayDataSource(vector7, sort_order="none")
#        value_range.add(y7)
#        value_mapper = LinearMapper(range=value_range)
#        
#        self.line_plot = LinePlot(index = index, value = y7,
#                                        index_mapper = index_mapper,
#                                        value_mapper = value_mapper,
#                                        color = tuple(COLOR_PALETTE[i]),
#                                        width = 20,
#                                        edge_color = 'blue',
#                                        border_visible = False)
#    
#        add_default_grids(self.line_plot)
#        add_default_axes(self.line_plot, vtitle="K_a", htitle="x")
#    
#        self.plot_container.add(self.line_plot)

#*****************************************************************************************************************************


    #            pan = PanTool(line_plot)
    #            zoom = SimpleZoom(line_plot, tool_mode="box", always_on=False)
    #            broadcaster.tools.append(pan)
    #            broadcaster.tools.append(zoom)
    
        # Add the traits inspector tool to the container
        #
     #           self.plot_container.tools.append(TraitsTool( self.plot_container ))
        
        self.line_plot.tools.append(PanTool(self.line_plot))
        self.line_plot.overlays.append(ZoomTool(self.line_plot))
    
     #   pd.set_data("y" + str(i), ydata[:,i])
    
        
        # change the color of the curves
#        i = i +1  
#        # Legend
#        plots["Label %d"%i] = self.line_plot
#
#        legend = Legend(component=self.plot_container, padding=10, align="ul")
#        legend.tools.append(LegendTool(legend, drag_button="right"))
#        self.plot_container.overlays.append(legend)
#
#        # Set the list of plots on the legend
#        legend.plots = plots
#   
#
#        # Add the title at the top
#        self.plot_container.overlays.append(PlotLabel("FORCE vs DISPLACEMENT",
#                                  component=self.plot_container,
#                                  font = "swiss 16",
#                                  overlay_position="top"))

 #       self.splot = Plot(pd)
    #    plot_s = Plot(pd, padding=50)
    
#        self.splot.plot(("index", "y_name"), color="red")

    #---------------------------------------------------------------------------
    #  Creates the table editing tool bar:
    #---------------------------------------------------------------------------

    def _create_toolbar ( self, parent ):
        """ Creates the table editing toolbar.
        """
        factory = self.factory
            
        panel = wx.Panel( parent, -1, style=wx.CLIP_CHILDREN )
        toolbar = MFnChacoEditorToolbar( parent = parent, editor = self )
        tb_sizer = wx.BoxSizer( wx.HORIZONTAL )
        panel.SetSizer(tb_sizer)

        self.toolbar = toolbar
        tb_sizer.Add( toolbar.control, 0 )
        tb_sizer.Add( ( 1, 1 ), 1, wx.EXPAND )

        return panel
    
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

#                factory  = self.factory
#                plotitem = factory.plotitem
#                x_values = getattr(self.object, plotitem.index)
#                y_values = getattr(self.object, plotitem.value)
                x_values = self.value.xdata
                y_values = self.value.ydata
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

               # plot_gc = PlotGraphicsContext(self._plot.outer_bounds)
               # plot_gc.render_component(self._plot)
                              
                #self._plot_container.outer_bounds = list((800,600))
#                plot_gc = PlotGraphicsContext((400,300), dpi=72.0)
#                plot_gc.render_component(self._plot_container)

#                self.line_plot.bounds = [500,300]
#                self.line_plot.padding = 50
#
#                win_size = self.line_plot.outer_bounds
#             #   win_size = self.component.outer_bounds
#                plot_gc = PlotGraphicsContext(win_size)
#        
#        
#                # Have the plot component into it
#                plot_gc.render_component(self.line_plot)
#
#                # Finally, we tell the graphics context to save itself to disk as an image.
#                plot_gc.save(path)

                DPI = 70.0
                size=(550,450)
             #   self.plot_container = create_plot()
                self.plot_container.bounds = list(size)
                self.plot_container.do_layout(force=True)
                
                gc = PlotGraphicsContext(size, dpi=DPI)
               

                #gc = GraphicsContext((size[0]+1, size[1]+1))
                gc.render_component(self.plot_container)
                gc.save(path)

                
            except:
                print "Error saving!"
                raise
            print "Plot saved."
        return

#-------------------------------------------------------------------------------
#  'MFnChacoEditorToolbar' class:
#-------------------------------------------------------------------------------
class MFnChacoEditorToolbar ( HasPrivateTraits ):
    """ Toolbar displayed in table editors.
    """
    #---------------------------------------------------------------------------
    #  Trait definitions:
    #---------------------------------------------------------------------------

    # Do not sort columns:
    save_data = Instance( Action,
                        {'name':    'Save as data',
                         'tooltip': 'Save the function values',
                         'action':  'on_savedata',
                         'enabled': True,
                         'image':   ImageResource( 'add' ) } )

     # Move current object up one row:
    save_fig = Instance( Action,
                        { 'name':    'Save as fig',
                         'tooltip': 'Save as figure',
                         'action':  'on_savefig',
                         'enabled': True,
                         'image':   ImageResource( 'save' ) } )

    # The table editor that this is the toolbar for:
    editor = Instance( _MFnChacoEditor )

    # The toolbar control:
    control = Any

    #---------------------------------------------------------------------------
    #  Initializes the toolbar for a specified window:
    #---------------------------------------------------------------------------

    def __init__ ( self, parent = None, **traits ):
        super( MFnChacoEditorToolbar, self ).__init__( **traits )
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


class MFnChacoEditor( BasicEditorFactory ):
    """ Editor factory for plot editors.
    """
    
    klass = _MFnChacoEditor
    
    adapter = Instance( MFnPlotAdapter )
    def _adapter_default(self):
        return MFnPlotAdapter()
