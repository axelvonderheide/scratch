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
#gjkgggggggg
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
        x = [0.30, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.40, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.50, 0.51, 0.52, 0.53, 0.54, 0.55 ]
        index = ArrayDataSource(x)
        index_range = DataRange1D()
        index_range.add(index)
        index_mapper = LinearMapper(range=index_range)
            
        value_range = DataRange1D(  low_setting = 0.0  )
        
            
#        vector1 = [0.3684, 0.3789, 0.3894, 0.4000, 0.4105, 0.4210, 0.4316, 0.4421, 0.4526, 0.4631, 0.4737, 0.4842, 0.4947, 0.5052, 0.5158, 0.5263, 0.5368, 0.5473, 0.5579, 0.5684, 0.5789, 0.5894, 0.6000, 0.6105, 0.6210, 0.6316]
#        y1 = ArrayDataSource(vector1, sort_order="none")
#        value_range.add(y1)
#        value_mapper = LinearMapper(range=value_range)
#        
#        self.line_plot = LinePlot(index = index, value = y1,
#                                        index_mapper = index_mapper,
#                                        value_mapper = value_mapper,
#                                        color = tuple(COLOR_PALETTE[i]),
#                                        width = 20,
#                                        edge_color = 'blue',
#                                        border_visible = False)
#    
#        add_default_grids(self.line_plot)
#        add_default_axes(self.line_plot, vtitle="d2/d=0.05   w1", htitle="x")
#    
#        self.plot_container.add(self.line_plot)
        

#        omega = [0.0041, 0.0146, 0.0252, 0.0357, 0.0462, 0.0567, 0.0673, 0.0778, 0.0883, 0.0988, 0.1094, 0.1199, 0.1304, 0.1410, 0.1515, 0.1620, 0.1725, 0.1831, 0.1936, 0.2041, 0.2146, 0.2252, 0.2357, 0.2462, 0.2567, 0.2673]
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
#        add_default_axes(self.line_plot, vtitle="d2/d=0.05   w2", htitle="x")
#    
#        self.plot_container.add(self.line_plot)
        
        
#        vector2 = [0.3686, 0.3797, 0.3908, 0.4020, 0.4131, 0.4242, 0.4353, 0.4464, 0.4575, 0.4686, 0.4797, 0.4908, 0.5020, 0.5131, 0.5242, 0.5353, 0.5464, 0.5575, 0.5686, 0.5797, 0.5908, 0.6020, 0.6131, 0.6242, 0.6353, 0.6464]
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
#        add_default_axes(self.line_plot, vtitle="d2/d=0.10   w1", htitle="x")
#        self.plot_container.add(self.line_plot)


#        vector3 = [0.0043, 0.0154, 0.0266, 0.0377, 0.0488, 0.0599, 0.0710, 0.0821, 0.0932, 0.1043, 0.1154, 0.1266, 0.1377, 0.1488, 0.1599, 0.1710, 0.1821, 0.1932, 0.2043, 0.2154, 0.2266, 0.2377, 0.2488, 0.2599, 0.2710, 0.2821]
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
#        add_default_axes(self.line_plot, vtitle="d2/d=0.10   w2", htitle="x")
#    
#        self.plot_container.add(self.line_plot)


#        vector4 = [0.3689, 0.3806, 0.3924, 0.4042, 0.4159, 0.4277, 0.4395, 0.4512, 0.4630, 0.4748, 0.4865, 0.4983, 0.5101, 0.5218, 0.5336, 0.5453, 0.5571, 0.5689, 0.5806, 0.5924, 0.6042, 0.6159, 0.6277, 0.6395, 0.6512, 0.6630 ]
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
#        add_default_axes(self.line_plot, vtitle="d2/d=0.15   w1", htitle="x")
#    
#        self.plot_container.add(self.line_plot)

#        vector5 = [0.0046, 0.0164, 0.0281, 0.0399, 0.0517, 0.0634, 0.0752, 0.0869, 0.0987, 0.1105, 0.1222, 0.1340, 0.1458, 0.1575, 0.1693, 0.1811, 0.1928, 0.2046, 0.2164, 0.2281, 0.2399, 0.2517, 0.2634, 0.2752, 0.2869, 0.2987]
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
#        add_default_axes(self.line_plot, vtitle="d2/d=0.15   w2", htitle="x")
#        self.plot_container.add(self.line_plot)


        vector6 = [0.3692, 0.3817, 0.3942, 0.4067, 0.4192, 0.4317, 0.4442, 0.4567, 0.4692, 0.4817, 0.4942, 0.5067, 0.5192, 0.5317, 0.5442, 0.5567, 0.5692, 0.5817, 0.5942, 0.6067, 0.6192, 0.6317, 0.6442, 0.6567, 0.6692, 0.6817]
        y6 = ArrayDataSource(vector6, sort_order="none")
        value_range.add(y6)
        value_mapper = LinearMapper(range=value_range)
        
        self.line_plot = LinePlot(index = index, value = y6,
                                        index_mapper = index_mapper,
                                        value_mapper = value_mapper,
                                        color = tuple(COLOR_PALETTE[i]),
                                        width = 20,
                                        edge_color = 'blue',
                                        border_visible = False)
    
        add_default_grids(self.line_plot)
        add_default_axes(self.line_plot, vtitle="d2/d=0.20   w1", htitle="x")
    
        self.plot_container.add(self.line_plot)
#
#        vector7 = [0.0055, 0.0194, 0.0334, 0.0474, 0.0614, 0.0753, 0.0893, 0.1033, 0.1173, 0.1312, 0.1452, 0.1592, 0.1732, 0.1871, 0.2011, 0.2151, 0.2291, 0.2430, 0.2570, 0.2710, 0.2850, 0.2989, 0.3129, 0.3269, 0.3409, 0.3548]
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
#        add_default_axes(self.line_plot, vtitle="d2/d=0.20   w2", htitle="x")
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
