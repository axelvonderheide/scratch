'''
Created on Jun 8, 2011

@author: schmerl
'''

from enthought.mayavi.core.api import PipelineBase
from enthought.mayavi.core.ui.api import MayaviScene, SceneEditor, \
    MlabSceneModel

from Data_Modell import DataModell
from enthought.traits.api import HasTraits, Range, Instance, on_trait_change, \
    Trait, Property, Constant, DelegatesTo, Str
from enthought.traits.ui.api import View, Item, Group

from enthought.mayavi import mlab




class threeDview( HasTraits ):
    
    #
    data = Instance(DataModell)
    
    
    
    
    scene = Instance( MlabSceneModel, () )
    plot = Instance( PipelineBase )
    
    
    # When the scene is activated
    @on_trait_change( 'scene.activated' )
    def update_plot( self ):
        # Create the mesh
        x = self.data.x
        y = self.data.y
        z = self.data.z
        
        #Color
        mlab.figure(1, fgcolor=(0, 0, 0), bgcolor=(1, 1, 1))
        #create Points
        pts = mlab.points3d(x, y, z, z, scale_mode='none', scale_factor=0)
        #Create mesh with nearest neighbors
        #TODO BETTER
        # Mesh has to be generated with crease line information
        mesh = mlab.pipeline.delaunay2d(pts)
        
        
        if self.plot is None:
            #Plott steuern
            self.plot = self.scene.mlab.pipeline.surface(mesh)
            self.plot.module_manager.scalar_lut_manager.label_text_property.italic = False
            
            
    view = View( Item( 'scene', editor = SceneEditor( scene_class = MayaviScene ),
                height = 250, width = 300, show_label = False ),
            Group( 
                     
                     ),
            resizable = True,
            )


if __name__ == '__main__':

    
    my_model = threeDview(data = DataModell())
    my_model.configure_traits()
