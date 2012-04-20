
from enthought.traits.api import \
    HasTraits, Float, Int, Array, Property, cached_property, \
    Tuple, List, Str, on_trait_change, Button, Delegate, \
    Instance, Trait

from enthought.traits.ui.api import \
    View, Item, Group, HGroup, VGroup, VSplit, HSplit, CheckListEditor, TextEditor
        
from math import floor
from numpy import zeros, mgrid, c_, indices, transpose, array, arange, \
    asarray, ix_, ones, random

# tvtk related imports
#
from enthought.traits.ui.api import View, Item, HSplit, VSplit, InstanceEditor
from enthought.tvtk.api import tvtk

from enthought.tvtk.pyface.scene_editor import SceneEditor
from enthought.mayavi.core.source import Source

#---------------------------------------------------------------------
# C L A S S  GeoNDGrid
#---------------------------------------------------------------------
class QuadSource(Source):

    height = Float( 1 )
    width = Float( 1 )

    #------------------------------------------------------------------------
    # Implementation of the mayavi source interface
    #------------------------------------------------------------------------    
    poly_data = Instance( tvtk.PolyData )
    def _poly_data_default(self):
        return tvtk.PolyData()
    
    def _outputs_default(self):
        return [ self._get_data() ]
    
    changed = Button('Draw')
    @on_trait_change('changed')
    def redraw(self):
        '''
        Rebuild the poly data for the current values
        '''
        self.data_changed = True
        self.outputs = [self._get_data()]

    def _get_data(self):
        points   = self._get_points()
        lines    = self._get_lines()
        faces    = self._get_faces()

        self.poly_data = tvtk.PolyData(points = points, lines = lines, polys = faces)
        arr = self._get_scalars() 
        if len( arr ) > 0:
            self.poly_data.point_data.scalars = arr
            self.poly_data.point_data.scalars.name = 'weibull random field'
        return self.poly_data
    
    def _get_points(self):
        '''
        Reshape the grid into a column.
        '''
        return array([[0,0,0],
                      [self.width/2.,   0,             0],
                      [self.width/2.,   self.height/2., 0],
                      [0,              self.height/2., 0],
                      [self.width,     0,             0],
                      [self.width,     self.height/2., 0],
                      [self.width/2.,   self.height,   0],
                      [self.width,     self.height,0]], 'f')
        
    def _get_lines(self):
        '''
        Get the number of lines.
        '''
        return array([[0,1],
                      [1,2],
                      [2,3],
                      [3,0]])


    def _get_faces(self):
        '''
        Only return data of n_dims = 2.
        '''
        return array([[0,3,2,1],[1,4,5,2],[2,5,7,6]])

    def _get_scalars(self):
        return random.weibull( 1, size = 8 ) 
    
    traits_view = View( HSplit( Group(Item( 'changed', show_label = False ), 
                                      Item('height', resizable = False ),
                                      Item('width', resizable = False )
                                      ),
                                            ),                        
                        resizable = True )


if __name__ == '__main__':

    from enthought.mayavi.scripts import mayavi2
    @mayavi2.standalone
    def view():
        from enthought.mayavi.modules.api import Outline, Surface
        # 'mayavi' is always defined on the interpreter.
        mayavi.new_scene()
        # Make the data and add it to the pipeline.
        mfn = QuadSource()
            
        mayavi.add_source(mfn)
        # Visualize the data.
        mayavi.add_module(Outline())
        mayavi.add_module(Surface())    
        
        
    view()
