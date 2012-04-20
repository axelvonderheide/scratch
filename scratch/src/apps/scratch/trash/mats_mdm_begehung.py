
from enthought.traits.api import \
     Array, Bool, Callable, Enum, Float, HasTraits, \
     Instance, Int, Trait, Range, HasTraits, on_trait_change, Event, \
     implements, Dict, Property, cached_property, Str

# Enthought library imports.
from enthought.pyface.api import ApplicationWindow, GUI
from enthought.pyface.action.api import Action, MenuManager, MenuBarManager
from enthought.pyface.action.api import StatusBarManager, ToolBarManager
from enthought.pyface.api import ApplicationWindow, SplitPanel

from mats_mdm2d_phi import MPArrayCMDM
from mats_mdm2d import MACMDMExplore

import wx

class MainWindow(ApplicationWindow):
    """ The main application window. """

    #### 'SplitApplicationWindow' interface ###################################

    # The direction in which the panel is split.
    title = Str('Composite-Microplane-Damage-Model-Explorer')

    ###########################################################################
    # 'object' interface.
    ###########################################################################

    def __init__(self, **traits):
        """ Creates a new application window. """

        # Base class constructor.
        super(MainWindow, self).__init__(**traits)

        # Create the window's menu, tool and status bars.
        self._create_action_bars()

        return

    def _create_action_bars(self):
        # Create an action that exits the application.
        exit_action = Action( name='E&xit',  on_perform=self.close )
        test_action = Action( name='S&tart', on_perform=self.start_calculation )

        # Add a menu bar.
        self.menu_bar_manager = MenuBarManager(
            MenuManager(exit_action, name='&File')
        )

        # Add some tool bars.
        self.tool_bar_managers = [
            ToolBarManager(
            test_action,
            exit_action, 
            name='Tool Bar 1', show_tool_names=True
            )
            ]

        # Add a status bar.
        self.status_bar_manager = StatusBarManager()
        self.status_bar_manager.message = 'Example application window'

    def _create_contents( self, parent ):

        panel = wx.Panel(parent, -1, style=wx.CLIP_CHILDREN)
        sizer = wx.BoxSizer(wx.VERTICAL)
        panel.SetSizer(sizer)
        panel.SetAutoLayout(True)
        
        self._mpe = MACMDMExplore()

        ui_alpha = self._mpe.edit_traits( parent = panel,
                                          view = 'traits_view_alpha',
                                          kind = 'subpanel' )

        ui_rv = self._mpe.sim.edit_traits( parent = panel,
                                           #view = 'traits_view_begehung',
                                           kind = 'subpanel' )
        
        sizer.Add(ui_alpha.control, 0, wx.EXPAND)
        sizer.Add(ui_rv.control, 1, wx.EXPAND)
        
        # Resize the panel to fit the sizer's minimum size.
        sizer.Fit(panel)

        # Set up the timer and start it up
        timerId = wx.NewId()
        self.timer = wx.Timer(panel, timerId)
        # Register a callback with the timer event
        panel.Bind(wx.EVT_TIMER, self._mpe.sim.tloop.rtrace_mngr.timer_tick, id=timerId)

        self._mpe.sim.tloop.KMAX = 1000
        self._mpe.sim.tloop.RESETMAX = 20
        self._mpe.sim.tloop.rtrace_mngr.timer = self.timer
        
        return panel

    def start_calculation( self ):
        self._mpe.sim.tloop.teval()

if __name__ == "__main__":
    # Create the GUI (this does NOT start the GUI event loop).
    gui = GUI()

    # Create and open the main window.
    window = MainWindow()
    window.open()
    window.size = (900,600)
    # Start the GUI event loop!
    gui.start_event_loop()
