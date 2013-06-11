


from enthought.mayavi import mlab
from enthought.mayavi.core.ui.api import MlabSceneModel, SceneEditor
from enthought.mayavi.core.ui.mayavi_scene import MayaviScene
from enthought.mayavi.tools.mlab_scene_model import MlabSceneModel
from enthought.traits.api import HasTraits, Float, Property, cached_property, \
    Event, Array, Instance, File, Int, Directory, Button, Range, Enum, \
     on_trait_change
from enthought.traits.ui.api import Item, View, Group, Handler, HSplit, VSplit, \
    HGroup, HSplit, VGroup, Tabbed, Label, Controller, ModelView
from enthought.traits.ui.api import View, Item, FileEditor, DirectoryEditor, \
    HistoryEditor
from enthought.traits.ui.menu import OKButton, CancelButton, Action, Menu, \
    MenuBar
from enthought.tvtk.pyface.scene_editor import SceneEditor 
from matplotlib.figure import Figure
from numpy import sin, cos, linspace, pi, mean, zeros
from traits.editors.mpl_figure_editor import MPLFigureEditor

from numpy.random import random


class YarnRawData( HasTraits ):
    '''
        Read raw data files
    '''
    raw_data_directory = Directory( 'directory', entries=10 )
    
    traits_view = View( Item( 'raw_data_directory' ) )
    
class YMBDataGrid( HasTraits ):
    def __call__( self ):
        pass
    
class YMBSlider( HasTraits ):
    '''
        Slicing data array (2d NumPy array)
        Return data for statistics (1d NumPy array)
    '''
    data = Instance( YarnRawData )
    evaluated_variable = Enum( 'radius', 'bond quality', 'bond free length', 'slack' )
    cut_slider = Range( 0, 10, mode='slider', auto_set=False, enter_set=True, modified=True )
    filament_slider = Range( 0, 1000, mode='slider', auto_set=False, enter_set=True, modified=True )

    
    traits_view = View( Item( 'evaluated_variable', springy=True ),
                        Item( 'cut_slider' , springy=True ),
                        Item( 'filament_slider', springy=True ) )
    
class Statistics( HasTraits ):
    yarn_slider = Instance( YMBSlider )
    def get_mean( self ):
        return mean( self.yarn_slider() )

    
class YMBHist( HasTraits ):
    raw_data = Instance( YarnRawData, () )
    slider = Instance( YMBSlider, () )
    figure = Instance( Figure )
    # The scene model.
    scene = Instance( MlabSceneModel, () )
    button = Button( 'Redraw' )
    
    def _figure_default( self ):
        figure = Figure()
        figure.add_axes( [0.1, 0.1, 0.8, 0.8] )
        return figure

    data_changed = Event( True )
    @on_trait_change( 'slider.+modified' )
    def _redraw( self ):
        figure = self.figure
        axes = figure.axes[0]
        axes.clear()
        axes.plot( [0, 5], [0, 5] )
        self.data_changed = True
    
    def _button_fired( self ):
        self.redraw_scene( self.scene )
        
    def redraw_scene( self, scene ):
        # Notice how each mlab call points explicitely to the figure it
        # applies to.
        mlab.clf( figure=scene.mayavi_scene )
        x, y, z, s = random( ( 4, 100 ) )
        mlab.points3d( x, y, z, s, figure=scene.mayavi_scene )
    
    view = View( Group( Item( 'raw_data', style='custom',
                                       show_label=False ),
                 Item( 'slider', style='custom', show_label=False ),
                HSplit( 
                Group( Item( 'button', show_label=False ),
                       Item( name='scene',
                     editor=SceneEditor( scene_class=MayaviScene ),
                     show_label=False ) ),
                Group( Item( 'figure', style='custom',
                              editor=MPLFigureEditor(),
                              show_label=False ) ),
                ),
                ),
                id='yarn_structure_view',
                resizable=True,
                        ) 
   
    
if __name__ == '__main__':
    data = YMBHist()
    data.configure_traits()
