
# Standard library imports.
import os, sys

# Put the Enthought library on the Python path.
sys.path.append(os.path.abspath(r'..\..\..'))

# Enthought library imports.
from enthought.pyface.api import FileDialog
from enthought.pyface.api import OK
from enthought.pyface.api import confirm, error, information, warning, YES
from enthought.pyface.api import ApplicationWindow, GUI
from enthought.pyface.action.api import Action, MenuBarManager, MenuManager, Separator
from enthought.pyface.tvtk.api import DecoratedScene
from enthought.traits.api import Instance
from enthought.pyface.action.api import StatusBarManager

from FourthStep import X2COSX
import pickle
from copy import copy

class NewAction(Action):

    def __init__(self, window):
        self._window = window
        self.name = "N&ew.."

    def perform(self):
        
        self._window.scene = X2COSX()
        self._window.scene.a = 5.0
        self._window.scene.b = 3.0
        self._window.scene.c = 1.0
        self._window.scene.refresh()

class SaveAction(Action):

    def __init__(self, window):
        self._window = window
        self.name = "S&ave.."

    def perform(self):
        extns = ['*.pickle']
        dlg = FileDialog(parent=self._window.control, action='save as',
                wildcard='|'.join(extns), title="Save As")
        if dlg.open() == OK:
            file = open( dlg.path, 'w' )
            pickle.dump(self._window.scene, file)
            file.close()  

class OpenAction(Action):
    def __init__(self, window):
        self._window = window
        self.name = "O&pen.."

    def perform(self):
        extns = ['*.pickle']
        dlg = FileDialog(parent=self._window.control, action='open',
                wildcard='|'.join(extns), title="Open")
        if dlg.open() == OK:
            file = open( dlg.path, 'r' )
            self._window.scene = pickle.load(file)
            self. refresh()
            file.close()

class MainWindow(ApplicationWindow):
    """ The main application window. """

    ###########################################################################
    # 'object' interface.
    ###########################################################################

    scene = Instance(X2COSX)

    def _create_contents(self, parent):
        
        self.scene = X2COSX()
        self.ui = self.scene.edit_traits(parent = parent, kind = 'panel')
        
        return self.ui.control

    def x__init__(self, **traits):
        """ Creates a new application window. """

        # Base class constructor.
        super(MainWindow, self).__init__(**traits)

        # Add a menu bar.
        self.menu_bar_manager = MenuBarManager(
            MenuManager(
                Action(name='E&xit', on_perform=self._on_exit),
                Separator(),
                NewAction(self),
                OpenAction(self),
                SaveAction(self),
                name = '&File',
            )  
        )

        # Add a status bar.
        self.status_bar_manager = StatusBarManager()
        self.status_bar_manager.message = 'Beispiel fuer Filedialog'
        
        return

    ###########################################################################
    # Private interface.
    ###########################################################################

    def _on_exit(self):

        parent = self.control

        if confirm(parent, 'exit?') == YES:
            self.close()

        return
    

# Application entry point.
if __name__ == '__main__':
    # Create the GUI (this does NOT start the GUI event loop).
    gui = GUI()

    # Create and open the main window.
    window = MainWindow(size=(500,500))
    window.open()

    # Start the GUI event loop!
    gui.start_event_loop()

##### EOF #####################################################################
