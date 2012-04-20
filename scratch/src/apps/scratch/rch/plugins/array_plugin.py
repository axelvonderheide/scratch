#!/usr/bin/env python
""" The entry point for an Envisage application. """

# Standard library imports.
import sys
import os.path
import logging

# Enthought library imports.
#from enthought.mayavi.plugins.app import get_plugins, setup_logger
from enthought.mayavi.plugins.app import setup_logger
from enthought.traits.api import List
from enthought.envisage.api import Plugin, ServiceOffer
from enthought.envisage.ui.workbench.api import WorkbenchApplication
from enthought.pyface.workbench.api import Perspective, PerspectiveItem

logger = logging.getLogger()

###############################################################################
# `spirridPerspective` class.
###############################################################################
class SpirridPerspective(Perspective):
    """ An default perspective for the app. """

    # The perspective's name.
    name = 'spirrid'

    # Should this perspective be enabled or not?
    enabled = True

    # Should the editor area be shown in this perspective?
    show_editor_area = True

    # View IDs.
    SPIRRID_VIEW = 'spirrid.spirrid' 

    # The contents of the perspective.
    contents = [
        PerspectiveItem(id=SPIRRID_VIEW, position='left'),
    ]

###############################################################################
# `spirridPlugin` class.
###############################################################################
class SpirridPlugin(Plugin):

    # Extension points we contribute to.
    PERSPECTIVES = 'enthought.envisage.ui.workbench.perspectives'
    VIEWS             = 'enthought.envisage.ui.workbench.views'
    SERVICE_OFFERS = 'enthought.envisage.ui.workbench.service_offers'

    # The plugin's unique identifier.
    id = 'spirrid.spirrid'

    # The plugin's name (suitable for displaying to the user).
    name = 'spirrid'

    # Perspectives.
    perspectives = List(contributes_to=PERSPECTIVES)
    
    # Services we contribute.
    service_offers = List(contributes_to=SERVICE_OFFERS)
    
    # Views.
    views = List(contributes_to=VIEWS)


    ######################################################################
    # Private methods.
    def _perspectives_default(self):
        """ Trait initializer. """
        return [SpirridPerspective]

    def _service_offers_default(self):
        """ Trait initializer. """
        spirrid_service_offer = ServiceOffer(
            protocol = 'array_container.SPIRRID',
            factory  = 'array_container.SPIRRID'
        )

        return [spirrid_service_offer]

    def _views_default(self):
        """ Trait initializer. """
        return [self._spirrid_view_factory]

    def _spirrid_view_factory(self, window, **traits):
        """ Factory method for spirrid views. """
        from enthought.pyface.workbench.traits_ui_view import \
                TraitsUIView

        spirrid = self._get_spirrid(window)
        tui_engine_view = TraitsUIView(obj=spirrid,
                                       id='spirrid.spirrid',
                                       name='spirrid',
                                       window=window,
                                       position='left',
                                       **traits
                                       )
        return tui_engine_view

    def _get_spirrid(self, window):
        """Return the spirrid service."""
        return window.get_service('array_container.SPIRRID')

def get_plugins():
    """Get list of default plugins to use for Mayavi."""
    from enthought.envisage.core_plugin import CorePlugin
    from enthought.envisage.ui.workbench.workbench_plugin import WorkbenchPlugin
    from enthought.plugins.python_shell.python_shell_plugin import PythonShellPlugin
    from enthought.plugins.text_editor.text_editor_plugin import TextEditorPlugin
    from enthought.tvtk.plugins.scene.scene_plugin import ScenePlugin
    from enthought.tvtk.plugins.scene.ui.scene_ui_plugin import SceneUIPlugin
    from enthought.mayavi.plugins.mayavi_plugin import MayaviPlugin
    from enthought.mayavi.plugins.mayavi_ui_plugin import MayaviUIPlugin
    from enthought.units.plugin.units_plugin import UnitsPlugin
#    from enthought.units.plugin.units_ui_plugin import UnitsUIPlugin
    
    plugins = [CorePlugin(),
               WorkbenchPlugin(),
               MayaviPlugin(),
               MayaviUIPlugin(),
               ScenePlugin(),
               SceneUIPlugin(),
               PythonShellPlugin(),
#              TextEditorPlugin()
            ]

    return plugins
######################################################################
def main():

    # Get the default mayavi plugins.
    plugins = get_plugins()

    # Inject our plugin up front so our perspective becomes the default.
    #plugins = [ SpirridPlugin() ]
    plugins.insert(0, SpirridPlugin())

    # Create an Envisage application.
    id = 'spirrid.spirrid'
    application = WorkbenchApplication(id=id,
                                       plugins = plugins 
                                       )
    # This needs to be done here since the ETSConfig.application_home is
    # not set correctly up to this point.
    setup_logger(logger, 'spirrid.log', mode=logging.ERROR)

    # Start the application.
    application.run()

# Application entry point.
if __name__ == '__main__':
    main()