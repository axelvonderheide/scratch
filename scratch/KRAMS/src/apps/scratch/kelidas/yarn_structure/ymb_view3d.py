
#import wxversion

#wxversion.select( '2.8' )


from enthought.mayavi import mlab
from enthought.mayavi.core.ui.api import MlabSceneModel, SceneEditor
from enthought.mayavi.core.ui.mayavi_scene import MayaviScene
from enthought.mayavi.tools.mlab_scene_model import MlabSceneModel
from enthought.traits.api import HasTraits, Float, Property, cached_property, \
    Event, Array, Instance, File, Int, Directory, Button, Range, Enum, \
     on_trait_change, Bool, Trait, HasPrivateTraits
from enthought.traits.ui.api import Item, View, Group, Handler, HSplit, VSplit, \
    HGroup, HSplit, VGroup, Tabbed, Label, Controller, ModelView
from enthought.traits.ui.api import View, Item, FileEditor, DirectoryEditor, \
    HistoryEditor
from enthought.traits.ui.menu import OKButton, CancelButton, Action, Menu, \
    MenuBar
from numpy import arange, pi, cos, sin, ones_like, array

from enthought.traits.api import HasTraits, Range, Instance, \
        on_trait_change, Trait
from enthought.traits.ui.api import View, Item, Group

from enthought.mayavi.core.api import PipelineBase
from enthought.mayavi.core.ui.api import MayaviScene, SceneEditor, \
                MlabSceneModel
from enthought.mayavi.tools.sources import MGlyphSource

from os.path import join
import os, re

from enthought.tvtk.api import tvtk
from enthought.mayavi.scripts import mayavi2

from enthought.mayavi import mlab
from ymb_data import YMBData

class YMBView3D( HasTraits ):
    # TODO : somewhere is small mistake MAG (1050,1080), probabily caused by zeros in connectivity data
    data = Instance( YMBData, () )

    var_enum = Trait( 'radius',
                      {'radius' : 'radius',
                       'shortest distance from the edge' : 'edge_distance',
                       'contact fraction' : 'cf',
                       'slack' : 'slack'} )
    scalar_arr = Property( depends_on='var_enum' )
    def _get_scalar_arr( self ):
        return getattr( self.data, self.var_enum_ )

    n_fiber = Property( depends_on='data' )
    start_fib = Range( 0, 2000, 0, mode='slider' )
    end_fib = Range( 0, 2000, 2000, mode='slider' )

    @cached_property
    def _get_n_fiber( self ):
        return 1000

    scene = Instance( MlabSceneModel, () )

    plot = Instance( PipelineBase )


    # When the scene is activated, or when the parameters are changed, we
    # update the plot.
    @on_trait_change( 'scalar_arr, start_fib, end_fib, scene.activated' )
    def update_plot( self ):
        x_arrr, y_arrr, z_arrr, r_arrr, d_arrr, cf_arrr = self.data.cut_data
        scalar_arrr = self.scalar_arr

        x_arr = x_arrr[self.start_fib:self.end_fib, :]
        y_arr = y_arrr[self.start_fib:self.end_fib, :]
        z_arr = z_arrr[self.start_fib:self.end_fib, :]
        scalar_arr = scalar_arrr[self.start_fib:self.end_fib, :]

        mask = y_arr > -1

        x = x_arr[mask]
        y = y_arr[mask]
        z = z_arr[mask]
        scalar = scalar_arr[mask]

        connections = -ones_like( x_arr )
        connections[mask] = range( 0, len( connections[mask] ) )
        connection = connections.copy()
        connection = connection.tolist()

        # TODO: better
        for i in range( 0, self.data.n_cuts + 1 ):
            for item in connection:
                try:
                    item.remove( -1 )
                except:
                    pass


        if self.plot is None:
            pts = self.scene.mlab.pipeline.scalar_scatter( x, y, z, scalar )
            pts.mlab_source.dataset.lines = connection
            self.plot = self.scene.mlab.pipeline.surface( 
                    self.scene.mlab.pipeline.tube( 
#                        self.scene.mlab.pipeline.stripper( 
                            pts,
#                        ),
                        tube_sides=10, tube_radius=15,
                    ),
                )

            self.plot.actor.mapper.interpolate_scalars_before_mapping = True
            self.plot.module_manager.scalar_lut_manager.show_scalar_bar = True
            self.plot.module_manager.scalar_lut_manager.show_legend = True
            self.plot.module_manager.scalar_lut_manager.shadow = True
            self.plot.module_manager.scalar_lut_manager.label_text_property.italic = False

        else:
            self.plot.mlab_source.dataset.reset()
            self.plot.mlab_source.dataset.points = array( [x, y, z] ).T
            self.plot.mlab_source.scalars = scalar
            self.plot.mlab_source.dataset.lines = connection




#        x, y, z, t = curve( self.n_meridional, self.n_longitudinal )
#        if self.plot is None:
#            self.plot = self.scene.mlab.plot3d( x, y, z, t,
#                                tube_radius=0.025, colormap='Spectral' )
#        else:
#            self.plot.mlab_source.set( x=x, y=y, z=z, scalars=t )


    # The layout of the dialog created
    view = View( Item( 'scene', editor=SceneEditor( scene_class=MayaviScene ),
                     height=250, width=300, show_label=False ),
                Group( 
                        '_', 'start_fib', 'end_fib', 'var_enum',
                     ),
                resizable=True,
                )




if __name__ == '__main__':
    my_model = YMBView3D()
    my_model.configure_traits()




