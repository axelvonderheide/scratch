
from enthought.traits.api import \
    HasTraits, Float, Int, Array, Property, cached_property, \
    Tuple, List, Str, on_trait_change, Button

from enthought.traits.ui.api import \
    View, Item, Group, HGroup, VGroup, VSplit, HSplit, CheckListEditor
    
# todo: MeshActors deleted - adapt for MVMeshSource    
from ibvpy.mesh.mv_mesh_source import MVMeshSource

from math import floor
from numpy import zeros, mgrid, c_, indices, transpose, array, arange, ix_, ones, random

# tvtk related imports
#
from enthought.traits.ui.api import View, Item, HSplit, VSplit, InstanceEditor
from enthought.tvtk.api import tvtk

from enthought.pyface.tvtk.scene_editor import SceneEditor
#---------------------------------------------------------------------
# C L A S S  MeshGrid
#---------------------------------------------------------------------
class MeshGrid(MVMeshSource):
    label = Str('MeshGrid')
    # map of coordinate labels to the indices
    _dim_map = {'x' : 0, 'y' : 1, 'z' :2 }
    
    # currently active dimensions
    active_dims = List( Str, ['x','y'] )

    # indices of the currently active dimensions
    dim_indices = Property( Array(int), depends_on = 'active_dims' )
    @cached_property
    def _get_dim_indices(self):
        ''' Get active indices '''
        return array( [ self._dim_map[dim_ix] for dim_ix in self.active_dims ], dtype='int_')

    # number of currently active dimensions
    n_dims = Property( Int, depends_on = 'active_dims' )
    @cached_property
    def _get_n_dims(self):
        return len( self.active_dims )

    # number of elements in each direction
    shape = List( Int, [1,1,1] )
    
    n_act_elems = Property( Array, depends_on = 'shape,dim_indices' )
    @cached_property
    def _get_n_act_elems(self):
        act_idx = ones( (3, ), int )
        shape = array( list(self.shape), dtype = int )
        act_idx[ self.dim_indices ] += shape[ self.dim_indices ]
        return act_idx
    
    # total number of nodes of the system grid
    n_nodes = Property( Int , depends_on = 'n_act_elems' )
    @cached_property
    def _get_n_nodes( self ):
        # number of nodes used for the geometry approximation
        return reduce( lambda x,y:x*y, self.n_act_elems )

    ### Get the array of the node numbers corresponding to the regular nodal grid:
    # ('_geo'): nodes of the system used for geometry approximation 
    enum_nodes = Property( Array, depends_on = 'n_nodes' )
    @cached_property
    def _get_enum_nodes(self):
        # Returns an array of element numbers respecting the grid structure
        # (the nodes are numbered first in x-direction, then in y-direction and
        # last in z-direction)
        return arange( self.n_nodes ).reshape( tuple( self.n_act_elems ) )

    grid = Property( Array, depends_on = 'n_dims,shape,dim_indices,x_mins,x_maxs' )
    @cached_property
    def _get_grid(self):
        # slice(start,stop,step) with step of type 'complex' leads to that number of divisions 
        # in that direction including 'stop' (see numpy: 'mgrid')
        slices = [ slice(x_min,x_max,complex(0,n_n))
                   for x_min,x_max,n_n in zip( self.x_mins, self.x_maxs, self.n_act_elems ) ]
        return mgrid[tuple( slices )] 

    # domain range given as a box with min and max coordinates
    x_mins = Array( Float, shape = (3,), value = [0.,0.,0.] )
    x_maxs = Array( Float, shape = (3,), value = [1.,1.,1 ] )
    
    changed = Button('Draw')
        
    def _get_points(self):
        # find out which dimensions are active and augment the inactive with zeros
        return c_[ tuple( [ self.grid[i].flatten() for i in range(3) ] ) ]

    def _get_n_lines(self):
        act_idx = ones( (3, ), int )
        act_idx[ self.dim_indices ] += self.shape[ self.dim_indices ]
        return reduce( lambda x,y: x*y, act_idx ) 
    
    def _get_lines(self):
        '''
        Only return data if n_dims = 1
        '''
        if self.n_dims != 1: return array([],int)

        en = self.enum_nodes
        
        tidx = ones( (3,), dtype='int_' )
        tidx[ self.dim_indices ] = -1
        slices = tuple( [ slice(0,idx) for idx in tidx ] )

        offsets = array( [en[0,0,0],en[0,1,0]], dtype = 'int_')
            
        base_node_list = self.enum_nodes[slices].flatten()
        lines = []
        for n in base_node_list:
            lines.append( list( offsets + n) )
        print lines
        return lines
        
    def _get_n_faces(self):
        act_idx = ones( (3, ), int )
        shape = array( self.shape, dtype = int )
        act_idx[ self.dim_indices ] = 0.0
        act_idx[ self.dim_indices ] += shape[ self.dim_indices ]
        return reduce( lambda x,y: x*y, act_idx ) 
    
    def _get_faces(self):
        '''
        Only return data of n_dims = 2.
        '''
        if self.n_dims != 2: return array([],int)
        
        # get the slices extracting all corner nodes with the smallest 
        # node number within the element
        # 
        tidx = ones( (3,), dtype='int_' )
        tidx[ self.dim_indices ] = -1
        slices = tuple( [ slice(0,idx) for idx in tidx ] )
        base_node_list = self.enum_nodes[slices].flatten()

        ijk_arr = zeros( (3,4), dtype = int )
        ijk_arr[ self.dim_indices[0] ] = [0,0,1,1]
        ijk_arr[ self.dim_indices[1] ] = [0,1,1,0]

        offsets = self.enum_nodes[ ijk_arr[0], ijk_arr[1], ijk_arr[2] ]

        n_faces = self._get_n_faces()
        faces = zeros( (n_faces, 4), dtype = 'int_' )

        for n, base_node in enumerate( base_node_list ):
            faces[n,:] = offsets + base_node
        return faces

    def _get_volumes(self):
        '''
        Only return data if ndims = 3
        '''
        if self.n_dims != 3: return array([],int)
        en = self.enum_nodes
        
        tidx = ones( (3,), dtype='int_' )
        tidx[ self.dim_indices ] = -1
        slices = tuple( [ slice(0,idx) for idx in tidx ] )

        offsets = array( [en[0,0,0],en[0,1,0],en[1,1,0],en[1,0,0],
                          en[0,0,1],en[0,1,1],en[1,1,1],en[1,0,1]], dtype = 'int_')
        base_node_list = self.enum_nodes[slices].flatten()
        
        n_faces = self._get_n_faces()
        faces = zeros( (n_faces, 8), dtype = 'int_' )
        for n in base_node_list: faces[n,:] = offsets + n

        return faces

    # needed
    var = Str('dummy')
    idx = Int(0)

    def _get_scalars(self):
        return random.weibull( 1, size = self.n_nodes ) 
    
    field_entity_type = 'quad'

    traits_view = View( HSplit( Group(Item( 'changed', show_label = False ), 
                                      Item( 'active_dims@', editor = CheckListEditor( values = [ 'x', 'y', 'z'], 
                                                                 cols   = 3 ) ),
                                      Item('x_mins@', resizable = True ),
                                      Item('x_maxs@' ),
                                      Item('shape@' ),
                                      ),
                               HSplit(
                                     VSplit(Item(name='engine_view',
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
                                     Item(name='scene', 
                                           editor=SceneEditor(),
                                           show_label=False,
                                           resizable=True,
                                           height=500,
                                           width=500),
                                           ),
                                            ),                        
                        resizable = True )

if __name__ == '__main__':
    mfn = MeshGrid( active_dims = ['x','y'], shape = [2,1,5], x_maxs = [5,5,15])
    mfn.configure_traits()
