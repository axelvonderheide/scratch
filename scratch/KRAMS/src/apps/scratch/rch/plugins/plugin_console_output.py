""" How to get scrollbar in a plugin view. """

# Enthought library imports.
from enthought.traits.api import Float, HasTraits, List, Instance, Button
from enthought.traits.ui.api import Item, HGroup, VGroup, View, Tabbed, Group, HSplit

import logging
# GLOBALS
logger = logging.getLogger()
logger.addHandler(logging.StreamHandler())
logger.setLevel(logging.ERROR)

# The dummy model defines just data traits.
#
class MyModel(HasTraits):
    """ The data view. """

    # The model that we are a view of.

    throw_exc = Button('Throw exception')
    def _throw_exc_fired(self):
        print 'XXXXXXXXXXXXXXXXXXXXXXXX'
        raise ValueError, 'This is exceptional'
    
    # The view traits.
    data1 = Float(1.0)
    data2 = Float(1.0)
    data3 = Float(1.0)
    data4 = Float(1.0)
    data5 = Float(1.0)
    data6 = Float(1.0)
    data7 = Float(1.0)
    data8 = Float(1.0)
    data9 = Float(1.0)
    data10= Float(1.0)
    data11= Float(1.0)
    data12= Float(1.0)
    data13= Float(1.0)
    data14= Float(1.0)
    data15= Float(1.0)
    data16= Float(1.0)
    data17= Float(1.0)
    data18= Float(1.0)
    data19= Float(1.0)
    data20= Float(1.0)

    traits_ui_view = View(
          HGroup( 
        Item('data1'),
        Item('data2'),
        Item('data3'),
        Item('data4'),
        Item('data5'),
        Item('data6'),
        Item('data7'),
        Item('data8'),
        Item('data9'),
        Item('data10'),
          ),
        Item('data11'),
        Item('data12'),
        Item('data13'),
        Item('data14'),
        Item('data15'),
        Item('data16'),
        Item('data17'),
        Item('data18'),
        Item('data19'),
        Item('data20'),
        Item('throw_exc'),
        id='mymodel.data',
        resizable=True,
        scrollable=True,
        height = 0.5,
        width = 0.6,
        dock = 'tab',
    )

# Prepare the perspective for my model
#
from enthought.pyface.workbench.api import Perspective, PerspectiveItem

class MyPerspective(Perspective):
    """ A perspective containing the default Lorenz views. """
    
    name             = 'My'
    show_editor_area = False

    contents = [
        PerspectiveItem(id='mymodel.data'),
    ]

# Construct a simple plugin for my model
#
from enthought.envisage.api import Plugin
from enthought.pyface.workbench.api import TraitsUIView
from enthought.envisage.core_plugin import CorePlugin
from enthought.envisage.ui.workbench.workbench_plugin import WorkbenchPlugin

class MyModelUIPlugin(Plugin):
    """ The MyModel UI plugin.

    This plugin is part of the 'MyModel' example application.

    """

    # Extension points Ids.
    PERSPECTIVES   = 'enthought.envisage.ui.workbench.perspectives'
    VIEWS          = 'enthought.envisage.ui.workbench.views'

    #### 'IPlugin' interface ##################################################

    # The plugin's unique identifier.
    id = 'mymodel.ui'

    # The plugin's name (suitable for displaying to the user).
    name = 'MyModel UI'

    # Perspectives.
    perspectives = List(contributes_to=PERSPECTIVES)
    def _perspectives_default(self):
        """ Trait initializer. """
        return [MyPerspective]
    
    # Views.
    views = List(contributes_to=VIEWS)
    def _views_default(self):
        """ Trait initializer. """
        return [self._create_data_view]

    ###########################################################################
    # Private interface.
    ###########################################################################

    def _create_data_view(self, **traits):
        """ Factory method for the data view. """

        mymodel_view = TraitsUIView(
            id   = 'mymodel.data',
            name = 'Data',
            obj  = MyModel(),
            **traits
        )

        return mymodel_view

#################################################################
# Three different contexts for the mymodel's view
#
# 1) start using the workbench
# now - make the wondow smaller - no scrollbar appears
#
import enthought.mayavi.plugins.app as mayavi_app

def start_using_workbench_application():
    """ Run the application. """

    from enthought.envisage.ui.workbench.api import WorkbenchApplication
    # Enthought plugins.
    from enthought.envisage.core_plugin import CorePlugin
    from enthought.envisage.developer.developer_plugin import DeveloperPlugin
    from enthought.envisage.developer.ui.developer_ui_plugin import DeveloperUIPlugin
    from enthought.envisage.ui.workbench.workbench_plugin import WorkbenchPlugin
    from enthought.envisage.ui.single_project.project_plugin import ProjectPlugin
    from enthought.envisage.ui.workbench.api import WorkbenchApplication
    from enthought.plugins.python_shell.python_shell_plugin import PythonShellPlugin

    plugins = mayavi_app.get_plugins()
    
    # Create an application with the specified plugins.
    application = WorkbenchApplication(
        plugins=[
#            CorePlugin(), WorkbenchPlugin(), 
            MyModelUIPlugin(),
        ] + plugins
    )

    # Run it! This starts the application, starts the GUI event loop, and when
    # that terminates, stops the application.
    application.run()

    return

if __name__ == '__main__':
    
    start_using_workbench_application()

#### EOF ######################################################################