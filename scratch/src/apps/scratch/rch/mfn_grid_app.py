# This code simulates something the user would like to do.  In this
# case the code allows a user to create a 3D cube of data (a numpy
# array), specify an equation for the scalars and view it using the
# mayavi plugin.  The only "envisage bits" are the code that let one
# grab the running mayavi instance and script it.  The application
# trait is set by Envisage and we use the application to get hold of
# the mayavi engine.  Then we show the data once the mayavi engine has
# started.

# Standard library imports.
import numpy
import scipy

# Enthought library imports
from enthought.traits.api import HasTraits, Button, Instance, \
     Any, Str, Array
from enthought.traits.ui.api import Item, View, TextEditor

from subsid.mathkit.mfn.mfn_grid.mfn_ndgrid import MFnNDGrid

######################################################################
# `Explorer3D` class.
######################################################################
class Explorer3D(HasTraits):
    """This class basically allows you to create a 3D cube of data (a
    numpy array), specify an equation for the scalars and view it
    using the mayavi plugin.
    """

    ########################################
    # Traits.
    
    # Set by envisage when this is offered as a service offer.
    window = Instance('enthought.pyface.workbench.api.WorkbenchWindow')

    function = Instance( MFnNDGrid )
    def _function_default(self):
        return MFnNDGrid(shape = (8,8,1), active_dims = ['x', 'y'] )                 

    ########################################
    # Our UI view.
    view = View(Item('function@'),
                resizable=True,
                scrollable=True,
                )

    def get_mayavi(self):
        from enthought.mayavi.plugins.script import Script
        return self.window.get_service(Script)

#    def _update_data_fired(self):
#        self.show_data()

    def _show_data(self):
        mayavi = self.get_mayavi()
        if mayavi.engine.current_scene is None:
            mayavi.new_scene()
        from enthought.mayavi.sources.array_source import ArraySource

        mayavi.add_source(self.function)

        mayavi = self.get_mayavi()
        from enthought.mayavi.modules.api import Outline, Surface        
        from enthought.mayavi.modules.axes import Axes
        # Visualize the data.
        o = Outline()
        mayavi.add_module(o)
        a = Axes()
        mayavi.add_module(a)
        s = Surface()
        mayavi.add_module(s)
        
    def _window_changed(self):
        m = self.get_mayavi()
        print 'setting data'
        if m.engine.running:
            self._show_data()
        else:
            # Show the data once the mayavi engine has started.
            m.engine.on_trait_change(self._show_data, 'started')
        
            