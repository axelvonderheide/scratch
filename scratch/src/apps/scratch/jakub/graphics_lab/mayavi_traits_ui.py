#!/usr/bin/env python

"""An example of how to create an almost complete mayavi UI inside a traits UI
view.  This does not use Envisage and provides a similar UI as seen in the full
mayavi application. 

"""

# Authors: Prabhu Ramachandran <prabhu [at] aero.iitb.ac.in>
# Copyright (c) 2007, Enthought, Inc.
# License: BSD Style.

# Standard imports.
from numpy import sqrt, sin, mgrid

# Enthought imports.
from enthought.traits.api import HasTraits, Instance, Property, Enum
from enthought.traits.ui.api import View, Item, HSplit, VSplit, InstanceEditor
from enthought.tvtk.pyface.scene_editor import SceneEditor 
from enthought.mayavi.view.engine_view import EngineView
from enthought.mayavi import mlab
from enthought.mayavi.tools.mlab_scene_model import MlabSceneModel
from enthought.mayavi.plugins.app import get_plugins
# Set mlab to use the simple backend instead of envisage.
mlab.options.backend = 'simple'

######################################################################
class Mayavi(HasTraits):

    # Generate data with mayavi or mlab.
    data_source = Enum('mayavi','mlab')

    # The scene model.
    scene = Instance(MlabSceneModel, ())

    # The mayavi engine view.
    engine_view = Instance(EngineView)

    # The current selection in the engine tree view.
    current_selection = Property


    ######################
    view = View(HSplit(VSplit(Item(name='engine_view',
                                   style='custom',
                                   resizable=True,
                                   show_label=False
                                   ),
                              Item(name='current_selection',
                                   editor=InstanceEditor(),
                                   enabled_when='current_selection is not None',
                                   style='custom', 
                                   springy=True,
                                   show_label=False),
                                   ),
                        VSplit(Item(name='data_source'),
                               Item(name='scene', 
                                    editor=SceneEditor(),
                                    show_label=False,
                                    resizable=True,
                                    height=500,
                                    width=500),
                               )
                        ),
                resizable=True,
                scrollable=True
                )

    def __init__(self, **traits):
        HasTraits.__init__(self, **traits)
        self.engine_view = EngineView(engine=self.scene.engine)

        # Hook up the current_selection to change when the one in the engine
        # changes.  This is probably unnecessary in Traits3 since you can show
        # the UI of a sub-object in T3.
        self.scene.engine.on_trait_change(self._selection_change,
                                          'current_selection')

        self._data_source_changed(self.data_source)

    def generate_data_mlab(self):
        """Demo's how to generate data from mayavi."""
        # Create some data
        X, Y = mgrid[-2:2:200j, -2:2:200j]
        R = 10*sqrt(X**2 + Y**2)
        Z = sin(R)/R

        # Here we are using mlab to generate data.  We could just as well have
        # used the mayavi API.
        self.scene.mlab.surf(X, Y, Z, colormap='gist_earth')

    def generate_data_mayavi(self):
        """Shows how you can generate data using mayavi instead of mlab."""
        from enthought.mayavi.sources.api import ParametricSurface
        from enthought.mayavi.modules.api import Outline, Surface 
        from enthought.mayavi.filters.api import WarpVector
        from enthought.mayavi.sources.vtk_data_source import VTKDataSource
        from enthought.tvtk.api import tvtk
        from numpy import array
        e = self.scene.engine
#        s = ParametricSurface()
#        e.add_source(s)
#        e.add_module(Outline())
#        e.add_module(Surface())
        # The numpy array data.
        #points = array([[0,0,0], [1,0,0], [0,1,0], [0,0,1]], 'f')
        points = array([[0,0,0], [1,0,0], [1,1,0], [0,1,0]], 'f')
        warp = array([[0,0,0], [100,0,0], [1,1,0], [0,1,0]])
        deformation = tvtk.DoubleArray()
        deformation.number_of_components = 3
        deformation.number_of_tuples = 4
        deformation.set_tuple3(0,0.,0.,0)
        deformation.set_tuple3(1,20.,-5.,0.)
        deformation.set_tuple3(2,15.,3.,0.)
        deformation.set_tuple3(3,-4.,2.,0)
        #triangles = array([[0,1,3], [0,3,2], [1,2,3], [0,2,1]])
        triangles = array([[0,1,2,3]])
        temperature = array([10., 20., -20., 10.])
        # The TVTK dataset.
        mesh = tvtk.PolyData(points=points, polys=triangles)
        #mesh = tvtk.UnstructuredGrid(points=points)
        #cel_type = 7
        #mesh.set_cells(cel_type, triangles)
        #mesh.point_data.scalars = temperature
        #mesh.point_data.scalars.name = 'Temperature'
        mesh.point_data.vectors = warp
        src = VTKDataSource(data = mesh)
        e.add_source(src)
        e.add_filter(WarpVector())
        e.add_module(Outline())
        e.add_module(Surface())
        
        
    def _selection_change(self, old, new):
        self.trait_property_changed('current_selection', old, new)

    def _get_current_selection(self):
        return self.scene.engine.current_selection

    def _data_source_changed(self, value):
        self.scene.mlab.clf()
        if value == 'mayavi':
            self.generate_data_mayavi()
        elif value == 'mlab':
            self.generate_data_mlab()


if __name__ == '__main__':
    get_plugins()
    m = Mayavi()
    m.configure_traits()
