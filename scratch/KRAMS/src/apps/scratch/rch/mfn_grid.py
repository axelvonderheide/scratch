
from enthought.traits.api import \
    HasTraits, Float, Int, Array, Property, cached_property, \
    Tuple, List, Str, on_trait_change, Button, Delegate, \
    Instance, Trait

from enthought.traits.ui.api import \
    View, Item, Group, HGroup, VGroup, VSplit, HSplit, CheckListEditor, TextEditor
        
from numpy import zeros, mgrid, c_, indices, transpose, array, arange, \
    asarray, ix_, ones, random, floor, outer, diagflat, vdot

# tvtk related imports
#
from enthought.traits.ui.api import View, Item, HSplit, VSplit, InstanceEditor
from enthought.tvtk.api import tvtk

from enthought.pyface.tvtk.scene_editor import SceneEditor
from enthought.mayavi.core.source import Source

from mesh_grid_source import MeshGrid, GridPoint

class MFnNDGrid( MeshGrid ):
    '''
    Grid extended with functional interface to set the point data explicitly.
    
    The point data can be set using the method set_values_in_range 
    and extracted using the get_value method for an arbitrary position 
    within the domain.  
    '''
    point_values = Array( float )
    def _point_values_default(self):
        return []
    
    @on_trait_change('shape,active_dims')
    def _reset_point_values(self):
        print 'resetting point values',self._get_n_act_nodes()
        self.point_values = zeros( self.n_act_nodes, dtype = float )
        
    def set_value_all( self, value ):
        ''' Set all entries to the specified value.'''
        self.point_values[:] = value

    def _get_idx(self, x):
        ''' Get index for a specified x,y coordinate '''

        n_act_elems = array( self.n_act_nodes[:], int )
        n_act_elems[ self.dim_indices ] -= 1
        idx = array( floor( ( array(x, float) / 
                             ( array(self.x_maxs[:], float) - 
                               array(self.x_mins[:], float) ) ) 
                               * n_act_elems ), int )
        return idx
    
    def set_value_in_range( self, value, x_min, x_max, min_max = 0):
        ''' Set entries within the specified box to the specified value '''
        idx_mins = self._get_idx( x_min )
        idx_mins[ self.dim_indices ] += 1
        idx_maxs = self._get_idx( x_max ) + 1
        slices = [ slice(idx_min,idx_max) 
                   for idx_min,idx_max in zip( idx_mins, idx_maxs ) ]
        self.point_values[ slices ] = value
        
    def get_value(self, x):
        '''Return interpolated value.
        '''
        # The interplation is done using bilinear shape functions. 
        # The numpy array functionality is exploited to locate 
        # the interpolation cell coordinates and the point values in 
        # its corner nodes.
        
        # get the idx of the bottom left corner 
        idx_mins = self._get_idx( x )
        # idx of the upper right corner
        idx_maxs = idx_mins + 1
        idx_maxs[ self.dim_indices ] += 1
        # construct the slice for the cell
        slices = array( idx_mins, dtype = 'object_' )
        slices[ self.dim_indices ] = array( [ slice(idx_mins[i], idx_maxs[i]) 
                                              for i in self.dim_indices ] )
        # get the corner cell values (note, this operation 
        # includes dimensional reduction, i.e. corresponding 
        # to the currently active indices
        cell_values = self.point_values[ tuple( slices ) ]

        # get slice for the coordinates of the lower left corner 
        lower_left_slices = ( tuple( self.dim_indices), ) + tuple(idx_mins)
        # delta of the specified x with respect to lower left corner
        dx = (array(x, float)[self.dim_indices] - self.grid[ lower_left_slices ])

        # get the slice of the upper right corner
        idx_mins1 = idx_mins
        idx_mins1[self.dim_indices] += 1
        upper_right_slices = ( tuple( self.dim_indices), ) + tuple(idx_mins1)
        lengths = self.grid[upper_right_slices] - self.grid[lower_left_slices]
        coeffs = dx / lengths

        # perform the interpolation using broadcasting
        shapes = diagflat( [1 for i in self.dim_indices] ) + 1
        c_list = [ array( [ 1 - coeffs[i], coeffs[i] ], float ).
                   reshape( tuple( shapes[i] ) ) for i in range(self.n_dims) ]
        N_fn = reduce( lambda x, y: x * y, c_list )
        value = vdot( N_fn, cell_values )
        return value

    def _get_scalars(self):
        return self.point_values.flatten()
    
from enthought.mayavi.modules.api import Outline, Surface
from enthought.mayavi.scripts import mayavi2
mayavi2.standalone(globals())
   
if __name__ == '__main__':
        
    # Tests to write - set a quadrangle and values in the corner nodes
    # check the returned values for all the corner nodes
    # 
    # Check the interpolation in four internal positions.

    # 'mayavi' is always defined on the interpreter.
    mayavi.new_scene()
    # Make the data and add it to the pipeline.
    mfn = MFnNDGrid( active_dims = ['x', 'y'], shape = (1, 1, 1), 
                    x_maxs = GridPoint( x = 5, y = 5, z = 15 ) )
    mfn.set_value_in_range( 4.3, [1,1,1], [5,5,2] )
    print 'value at (1.1,1.1,1.1)',mfn.get_value([1.1,1.1,1])
    mayavi.add_source(mfn)
    # Visualize the data.
    mayavi.add_module(Outline())
    mayavi.add_module(Surface())    
