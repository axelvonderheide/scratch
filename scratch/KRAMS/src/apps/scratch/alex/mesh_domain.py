from enthought.traits.api import \
     Array, Bool, Enum, Float, HasTraits, HasStrictTraits, \
     Instance, Int, Trait, Str, Enum, \
     Callable, List, TraitDict, Any, Range, \
     Delegate, Event, on_trait_change, Button, \
     Interface, WeakRef, implements, Property, cached_property, Tuple, \
     Dict
from enthought.traits.ui.api import Item, View, HGroup, ListEditor, VGroup, \
     HSplit, Group, Handler, VSplit, TableEditor, ListEditor

from enthought.traits.ui.menu import NoButtons, OKButton, CancelButton, \
     Action

from enthought.traits.ui.ui_editors.array_view_editor \
    import ArrayViewEditor

from enthought.traits.ui.table_column \
    import ObjectColumn, ExpressionColumn

from enthought.traits.ui.table_filter \
    import TableFilter, RuleTableFilter, RuleFilterTemplate, \
           MenuFilterTemplate, EvalFilterTemplate, EvalTableFilter

from numpy import ix_, mgrid, array, arange, c_, newaxis, setdiff1d

from core.rv import RTrace

# tvtk related imports
#
from enthought.traits.ui.api import View, Item, HSplit, VSplit, InstanceEditor
from enthought.tvtk.api import tvtk
from core.sctx import \
    ISDomain, SDomain

from enthought.pyface.tvtk.scene_editor import SceneEditor

def flatten( x ):
    """flatten(sequence) -> list

    Returns a single, flat list which contains all elements retrieved
    from the sequence and all recursively contained sub-sequences
    (iterables).

    Examples:
    >>> [1, 2, [3,4], (5,6)]
    [1, 2, [3, 4], (5, 6)]
    >>> flatten([[[1,2,3], (42,None)], [4,5], [6], 7, MyVector(8,9,10)])
    [1, 2, 3, 42, None, 4, 5, 6, 7, 8, 9, 10]"""

    result = []
    for el in x:
        #if isinstance(el, (list, tuple)):
        if hasattr( el, "__iter__" ) and not isinstance( el, basestring ):
            result.extend( flatten( el ) )
        else:
            result.append( el )
    return result


# The definition of the demo TableEditor:
node_list_editor = TableEditor( 
    columns = [ ObjectColumn( label = 'Id', name = 'id_number' ),
                ObjectColumn( label = 'DOFs', name = 'dofs' ),
                ],
    editable = False,
    )

#-------------------------------------------------------------------
# MENode - part of the spatial domain
#-------------------------------------------------------------------

class NodeDOF( HasTraits ):
    id_number = Int
    dofs = List
    # check if the node of the regular grid is used in the 
    # element formulation (cf. 'get_elements'). Use this node in the assignement of the 
    # dof-numbers only if the attribute 'killed' is 'False'
    # see (cf. 'get_n_dofs', 'get_sliced_dofs')
    killed = Bool( False )

# The definition of the demo TableEditor:
elem_list_editor = TableEditor( 
    columns = [ ObjectColumn( label = 'Id', name = 'id_number' ),
                ObjectColumn( label = 'nodes', name = 'nodes' ),
                ExpressionColumn( label = 'Coordinates',
                                  expression = "str( object.get_X_mtx() )" ),
                ExpressionColumn( label = 'DOFs',
                                  expression = "str( object.get_dof_map() )" )
                ],
    editable = False,
    )

#-------------------------------------------------------------------
# MElem - spatial domain of the finite element
#-------------------------------------------------------------------

class MElem( HasTraits ):
    '''
    Finite element spatial representation.
    '''
    id_number = Int
    nodes_geo = Array( Int )
    nodes_dof = Array( Int )
    active = Bool( True )
    domain = WeakRef( ISDomain )

    def get_X_mtx( self ):
        '''
        Index mapping from the global array of coordinates.
        '''
        return self.domain.X_mtx[ ix_( list( self.nodes_geo ) )]

    def get_dof_map( self ):
        '''
        Return the dof map for the current element as a list
        '''
        dofs = []
        for n in self.nodes_dof:
            dofs += self.domain.nodes_dof[n].dofs
        return dofs

    def get_node_map( self ):
        '''
        Return the node map for the current element as a list
        '''
        nodes = []
        for n in self.nodes_dof:
            nodes.append( self.domain.nodes_dof[n].id_number )
        return nodes

#-------------------------------------------------------------------
# MeshDomain - spatial domain for the discretized region
#-------------------------------------------------------------------

class MGridDomain( SDomain, MeshActors ):
    label = Str( 'MGridDomain' )

    implements( ISDomain )
    #-----------------------------------------------------------
    # Primary parameters
    #-----------------------------------------------------------


    # Lengths of the system in direction x,y,z, respectively
    lengths = Tuple( ( 1., 1., 1. ) )
#    def _lenghts_default(self):
#        print 'default lengths_default'
#        return Tuple(1.,1.,1. )


    # Number of elements of the system in direction x,y,z, respectively
    shape = Tuple( ( 1, 1, 1 ) )
#    def _shape_default(self):
#        print 'default shape_default'
#        return (1,1,1)


    # Number of element nodes is defined by the number of devisions in 
    # each coordinate direction,e.g. 'n_e_nodes_geo = Tuple(2,1,1)' leads
    # to a 3x2x2 nodal grid of the element different grids can be defined
    # for geometry and field approximation. 
    # NOTE: The operands of 'n_e_nodes_geo'must be greater then zero if 
    # elements are defined in this direction!
    n_e_nodes_geo = Tuple( ( 1, 1, 1 ) )
#    n_e_nodes_geo = Tuple
#    def _n_e_nodes_geo_default(self):
#        print 'n_e_nodes_geo_default'
#        return (1,1,1)

    n_e_nodes_dof = Tuple( ( 1, 1, 1 ) )
#    n_e_nodes_dof = Tuple
#    def _n_e_nodes_dof_default(self):
#        print 'n_e_nodes_dof_default'
#        return (1,1,1)

    node_map_dof = List( range( 0, 8 ) )
#    node_map_dof = List
#    def _node_map_dof_default(self):
#        print 'node_map_dof_default'
#        return range(0,8)

    node_map_geo = List( range( 0, 8 ) )
#    node_map_geo = List
#    def _node_map_geo(self):
#        print 'node_map_geo_default'
#        return range(0,8)

    n_nodal_dofs = Int( 1 )

    # list of the nodemumbers of the regular grid that are 
    # not used in the element formulation, i.e. for which NodeDof.killed = True

# @todo: (alex) bring this into the Verify-format from Traits
#    # check if the specified parameters do not contradict with each other:
#    check_param_spec = Property( Any, depends_on = 'shape, length, n_e_nodes_geo, n_e_nodes_dof' )
#    @cached_property
#    def _get_check_param_spec( self ):
#        shape = self.shape
#        n_e_nodes_dof = self.n_e_nodes_dof
#        n_e_nodes_geo = self.n_e_nodes_geo
#        lengths = self.lengths
#        # Check if the operands in 'n_e_nodes_geo' differ from zero in case 
#        # elements are defined in the corresponding  direction!
#        warning_msg = ""
#        for i in range(3):
#            # this leads to reduction to one node only in direction i
#            if (shape[i] != 0 and n_e_nodes_geo[i] == 0): 
#                warning_msg += '!!! WARNING !!!: n_e_nodes_geo['+str(i)+']=0 but shape['+str(i)+'] is NOT 0 ! \n'
#            elif (shape[i] != 0 and n_e_nodes_dof[i] == 0): 
#                warning_msg += '!!! WARNING !!!: n_e_nodes_dof['+str(i)+']=0 but shape['+str(i)+'] is NOT 0 ! \n'
#            # this leads to two nodes placed at the same position    
#            elif (shape[i] != 0 and lengths[i] == 0): 
#                warning_msg += '!!! WARNING !!!: lengths['+str(i)+']=0 but shape['+str(i)+'] is NOT 0 ! \n'
#            # this does not effect the calculation but check input for consitstency   
#            elif shape[i] == 0 and (lengths[i] != 0 or n_e_nodes[i] != 0): 
#                warning_msg += '!!! NOTE !!!: lengths['+str(i)+']!=0 or n_e_nodes['+str(i)+']!=0 but shape['+str(i)+'] is 0 ! Check for consistency!\n'
#        if len(warning_msg) != 0:
#            print warning_msg               
#        return


    #-----------------------------------------------------------
    # Derived attributes - mesh representation
    #-----------------------------------------------------------  
    # Example of using the cached property to handle the mesh generation

    # total number of nodes of the system (used for the geometry approximation '_geo')
    n_nodes_geo = Property( Int , depends_on = 'shape,n_e_nodes_geo' )
#    @cached_property
    def _get_n_nodes_geo( self ):
        # number of nodes used for the geometry approximation
        n_nodes_list_geo = [ self.shape[i] * self.n_e_nodes_geo[i] + 1 for i in range( 0, 3 ) ]
        n_nodes_geo = n_nodes_list_geo[0] * n_nodes_list_geo[1] * n_nodes_list_geo[2]
        #print 'get total number of nodes_geo', n_nodes_geo
        return n_nodes_geo

    # total number of nodes of the system (used for the field approximation '_dof')
    n_nodes_dof = Property( Int, depends_on = 'shape,n_e_nodes_dof' )
    @cached_property
    def _get_n_nodes_dof( self ):
        # number of nodes used for field approximation
        n_nodes_list_dof = [ self.shape[i] * self.n_e_nodes_dof[i] + 1 for i in range( 0, 3 ) ]
        n_nodes_dof = n_nodes_list_dof[0] * n_nodes_list_dof[1] * n_nodes_list_dof[2]
        print 'get total number of nodes_dof', n_nodes_dof
        return n_nodes_dof

    X_mtx = Property( Array, depends_on = 'shape,lengths,n_e_nodes_geo' )
    @cached_property
    def _get_X_mtx( self ):
        print 'setting up matrix of global coordinates'
        n_e_geo = [ self.shape[i] * self.n_e_nodes_geo[i] for i in range( 0, 3 ) ]
        # slice(start,stop,step) with step of type 'complex' leads to that number of divisions 
        # in that direction including 'stop' (see numpy: 'mgrid')
        slices = [ slice( 0, length, complex( 0, n_e + 1 ) )
                   # the call 'lengths[::-1]' changes the order of the operands in 'length'
                   for length, n_e in zip( self.lengths[::-1], n_e_geo[::-1] ) ]
        c = mgrid[tuple( slices )]
        z, y, x = c
        # @TODO - this always copies the array
        #
        # create an array out of the three indicated lists interpreted as columns 
        # (see numpy: 'c_[index_expression])'
        X_mtx = c_[x.flatten(), y.flatten(), z.flatten()]
        print 'MGridDomain.X_mtx', X_mtx
        return X_mtx

    def get_X_mtx_dof( self ):
        '''setting up matrix of global coordinates for the dof-nodes based 
           on the regular rectangular grid (corresponding to the list nodes_dof'''
        n_e_dof = [ self.shape[i] * self.n_e_nodes_dof[i] for i in range( 0, 3 ) ]
        # slice(start,stop,step) with step of type 'complex' leads to that number of divisions 
        # in that direction including 'stop' (see numpy: 'mgrid')
        slices = [ slice( 0, length, complex( 0, n_e + 1 ) )
                   # the call 'lengths[::-1]' changes the order of the operands in 'length'
                   for length, n_e in zip( self.lengths[::-1], n_e_dof[::-1] ) ]
        c = mgrid[tuple( slices )]
        z, y, x = c
        # @TODO - this always copies the array
        #
        # create an array out of the three indicated lists interpreted as columns 
        # (see numpy: 'c_[index_expression])'
        X_mtx_dof = c_[x.flatten(), y.flatten(), z.flatten()]
        return X_mtx_dof



    ### Get the array of the node numbers corresponding to the regular nodal grid:
    # ('_geo'): nodes of the system used for geometry approximation 
    enum_nodes_geo = Property( Array, depends_on = 'shape,n_e_nodes_geo' )
    @cached_property
    def _get_enum_nodes_geo( self ):
        print 'enumerate nodes_geo'
        n_nodes_geo = self.n_nodes_geo
        # Returns an array of element numbers respecting the grid structure
        # (the nodes are numbered first in x-direction, then in y-direction and
        # last in z-direction)
        n_nodes_list_geo = [ self.shape[i] * self.n_e_nodes_geo[i] + 1 for i in range( 0, 3 ) ]
        enumt = arange( n_nodes_geo ).reshape( n_nodes_list_geo[2],
                                              n_nodes_list_geo[1],
                                              n_nodes_list_geo[0] )
        return enumt.transpose()

    # ('_dof'): nodes of the system used for field approximation 
    enum_nodes_dof = Property( Array, depends_on = 'shape,n_e_nodes_dof' )
    @cached_property
    def _get_enum_nodes_dof( self ):
        print 'enum_nodes_dof'
        # Array of element numbers of the entire system respecting the grid structure
        #
        n_nodes_list_dof = [ self.shape[i] * self.n_e_nodes_dof[i] + 1 for i in range( 0, 3 ) ]
        enumt = arange( self.n_nodes_dof ).reshape( n_nodes_list_dof[2],
                                                   n_nodes_list_dof[1],
                                                   n_nodes_list_dof[0] )
        return enumt.transpose()


    nodes_dof = Property( List )
    @cached_property
    def _get_nodes_dof( self ):
        print 'setting up nodes_dof'
        return [ NodeDOF( id_number = n ) for n in range( self.n_nodes_dof ) ]


    nodes_geo = Property( List )
    @cached_property
    def _get_nodes_geo( self ):
        print 'setting up nodes_geo'
        return [ NodeDOF( id_number = n ) for n in range( self.n_nodes_geo ) ]

    shape = Property( Int, depends_on = 'shape' )
    @cached_property
    def _get_shape( self ):
        print 'shape - calculating the number of elements'
        return max( 1, self.shape[0] ) * max( 1, self.shape[1] ) * max( 1, self.shape[2] )

    # List of finite element 
    elements = Property( List, \
                         depends_on = 'n_e_nodes_geo, node_map_geo' + \
                                      'n_e_nodes_dof, node_map_dof, shape' )
    @cached_property
    def _get_elements( self ):
        '''
        This procedure emulates a trivial mesh generation on a regular basis. 
        It might be replaced with a different kind of mesh. All other methods 
        and visualization should remain.
        NOTE:
        A different nodal grid with a independent numbering of the nodes is applied for 
        the field approximation ('_dof') and the geometry approximation ('_geo'). 
        The method assigns the node numbers of the nodal grid to the corresponding elements 
        based on the order defined in the list 'node_map'
        '''
        print 'setting up elements'
        elist = []
        e = 0
        #
        # specify the limits of the loop in order to avoid an empty range
        range_x = range( max( 1, self.shape[0] ) )
        range_y = range( max( 1, self.shape[1] ) )
        range_z = range( max( 1, self.shape[2] ) )
        # loop over all elements:
        for i in range_x:
            for j in range_y:
                for k in range_z:

                    # '_geo': get the node numbers of the system corresponding to the element
                    nodes_geo_arr_t = self.enum_nodes_geo[self.n_e_nodes_geo[0] * i:self.n_e_nodes_geo[0] * ( i + 1 ) + 1,
                                                          self.n_e_nodes_geo[1] * j:self.n_e_nodes_geo[1] * ( j + 1 ) + 1,
                                                          self.n_e_nodes_geo[2] * k:self.n_e_nodes_geo[2] * ( k + 1 ) + 1]
                    nodes_geo_arr = nodes_geo_arr_t.transpose()
                    nodes_geo_flat = flatten( nodes_geo_arr )
                    nodes_geo = array( nodes_geo_flat )[ ix_( self.node_map_geo ) ]

                    # '_dof': get the node numbers of the system corresponding to the element
                    nodes_dof_arr_t = self.enum_nodes_dof[self.n_e_nodes_dof[0] * i:self.n_e_nodes_dof[0] * ( i + 1 ) + 1,
                                                          self.n_e_nodes_dof[1] * j:self.n_e_nodes_dof[1] * ( j + 1 ) + 1,
                                                          self.n_e_nodes_dof[2] * k:self.n_e_nodes_dof[2] * ( k + 1 ) + 1]
                    nodes_dof_arr = nodes_dof_arr_t.transpose()
                    nodes_dof_flat_grid = flatten( nodes_dof_arr )
                    nodes_dof = array( nodes_dof_flat_grid )[ ix_( self.node_map_dof ) ]
                    # assign the node numbers ('_dof' and '_geo') to the element                       
                    elist.append( MElem( id_number = e,
                                         nodes_geo = nodes_geo,
                                         nodes_dof = nodes_dof,
                                         domain = self ) )
                    e += 1
        return elist


    killed_nodes_list = Property( List, depends_on = 'shape, n_e_nodes_dof, node_map_dof' )
    @cached_property
    def _get_killed_nodes_list( self ):
        '''
        Collect all nodes that are not used in any element.
        '''
        # specify the limits of the loop in order to avoid an empty range
        range_x = range( max( 1, self.shape[0] ) )
        range_y = range( max( 1, self.shape[1] ) )
        range_z = range( max( 1, self.shape[2] ) )
        # check if some nodes in 'nodes_dof_list' are not assigned  
        # to the element by the defined node_map, in this case these nodes
        # are not used in the element formulation and the list 'node_dof_flat_grid'
        # contains some nodenumbers that are not contained in 'node_dof_flat' 
        if ( max( 1, self.shape[0] ) * self.n_e_nodes_dof[0] + 1 ) * \
           ( max( 1, self.shape[1] ) * self.n_e_nodes_dof[1] + 1 ) * \
           ( max( 1, self.shape[2] ) * self.n_e_nodes_dof[2] + 1 ) > len( self.node_map_dof ):
            killed_nodes_list = []
            # loop over all elements:
            for i in range_x:
                for j in range_y:
                    for k in range_z:
                        # '_dof': get the node numbers of the system corresponding to the element
                        nodes_dof_arr_t = self.enum_nodes_dof[self.n_e_nodes_dof[0] * i:self.n_e_nodes_dof[0] * ( i + 1 ) + 1,
                                                              self.n_e_nodes_dof[1] * j:self.n_e_nodes_dof[1] * ( j + 1 ) + 1,
                                                              self.n_e_nodes_dof[2] * k:self.n_e_nodes_dof[2] * ( k + 1 ) + 1]
                        nodes_dof_arr = nodes_dof_arr_t.transpose()
                        nodes_dof_flat_grid = flatten( nodes_dof_arr )
                        nodes_dof = array( nodes_dof_flat_grid )[ ix_( self.node_map_dof ) ]
                        nodes_dof_flat = flatten( nodes_dof )
                        # list of the nodes of the regular grid which are not used and therefore 
                        # need to be killed in order to prevent that they are used in the 
                        # enumeration of the dofs (i.e. in 'get_dofs') 
                        killed_nodes_array_e = setdiff1d( nodes_dof_flat_grid, nodes_dof_flat )
                        # add the nodenumbers of the killed nodes (of the current checked element) 
                        # to the global variable 'killed_nodes_list' 
                        killed_nodes_list += list( killed_nodes_array_e )
                    #print 'self.killed_nodes_list',killed_nodes_list
        return killed_nodes_list


    # Number of dofs is a cached property of the FELineDomain. It must be
    # recalculated whenever the elements list gets changed.
    #
    n_dofs = Property( Int, depends_on = 'n_nodal_dofs' )

    # Update n_dofs if the element list has changed.  The
    # counting of dofs includes their enumeration as well.
    #
    @cached_property
    def _get_n_dofs( self ):
        print 'enumerating DOFs'
        # mark the Nodes in the killed list as killed
        for knode in self.killed_nodes_list:
            self.nodes_dof[knode].killed = True
            print 'self.nodes_dof[knode].killed', self.nodes_dof[knode].killed

        n_dofs = 0
        for node in self.nodes_dof:
            # check if some of the nodes in 'nodes_dof' are to be killed
            # because they are not used in the element formulation (e.g. serendipity elements)
            if node.killed == False:
                # assign the dofs to the node in the order used in the element formulation
                # (i.e. the order of the shape functions)
                node.dofs = range( n_dofs, n_dofs + self.n_nodal_dofs )
                print 'node.dofs', node.dofs
                print 'node.killed', node.killed
                n_dofs += self.n_nodal_dofs
        print 'n_dofs', n_dofs
        return n_dofs


    def get_sliced_dofs( self, slice0, slice1, slice2 ):
        '''
        returns a two dimensional array of dof-numbers based on the selected nodes 
        of the regular grid. The rows correspond to the selected nodes and the 
        columns to dofs in x-,y-, and z-direction respectively
        e.g. for a linear brick element "get_left_dofs()" yields
        [[ 0  1  2]
         [12 13 14]
         [ 6  7  8]
         [18 19 20]]        
        To constrain the left dofs in x-direction slicing can be applied on the returned array, e.g
        get_left_dofs()[:,0]. To constraint the left dofs in z-direction the call is: 
        get_left_dofs()[:,2]. To constraint only a single dof slicing can also be applied, e.g 
        get_left_dofs()[3,1] constraints node three in nodes_list in y-direction.    
        '''
        # 'enum_nodes' and 'n_dofs' are cached properties
        # access to 'n_dofs' initiates the enumeration of the dofs
        # therefore, it must be invoked before accessing 
        # the dofs in the nodes.
        # @todo: 
        self.n_dofs
        #
        nodes_list_grid = flatten( self.enum_nodes_dof[ slice0, slice1, slice2 ] )
        #
        # Remove the killed nodes from the list. The killed nodes have no attribute 'dof'
        nodes_list = flatten( setdiff1d( nodes_list_grid, self.killed_nodes_list ) )
#        print 'nodes_list_grid', nodes_list_grid
#        print 'nodes_list', nodes_list
#        print 'self.killed_nodes_list', self.killed_nodes_list        
        #
        dofs_list = [self.nodes_dof[n].dofs for n in nodes_list]
        print 'dofs_list', dofs_list
        return array( dofs_list, dtype = 'int_' )

    ### Some predefined definitions using 'get_sliced_dofs':
    # use of 'slice(None)' coorsponds to 'slice(None,None,None)',
    # i.e. select all values (i.e. the nodenumbers) in this axis 
    #  0  - selects the first nodes of the grid 'enum_dof' in the 
    #       correspoding coordinate-direction
    #(-1) - selects the last nodes of the grid 'enum_dof' in the 
    #       correspoding coordinate-direction

    def get_left_dofs( self ):
       return self.get_sliced_dofs( 0, slice( None ), slice( None ) )

    def get_right_dofs( self ):
        return self.get_sliced_dofs( -1, slice( None ), slice( None ) )

    def get_bottom_dofs( self ):
        return self.get_sliced_dofs( slice( None ), 0, slice( None ) )

    def get_top_dofs( self ):
        return self.get_sliced_dofs( slice( None ), -1, slice( None ) )

    def get_back_dofs( self ):
        return self.get_sliced_dofs( slice( None ), slice( None ), 0 )

    def get_front_dofs( self ):
        return self.get_sliced_dofs( slice( None ), slice( None ), -1 )

    def get_top_middle_dofs( self ):
        if self.shape[0] % 2 == 0:
            # odd element number in x-direction (select only one node):
            slice_middle_x = self.shape[0] * self.n_e_nodes_dof[0] / 2
        else:
            print 'Error in get_top_middle_dofs: the method is only defined for an even number of elements in x-direction'
        return self.get_sliced_dofs( slice_middle_x, -1, slice( None ) )

    def get_bottom_left_dofs( self ):
        return self.get_sliced_dofs( 0, 0, slice( None ) )

    def get_bottom_right_dofs( self ):
        return self.get_sliced_dofs( -1, 0, slice( None ) )


    #---------------------------------------------------------------
    # Visualization related methods
    #---------------------------------------------------------------
    def _get_points( self ):
        return self.X_mtx
    def _get_point_numbers( self ):
        return arange( self.n_nodes_geo )
    def _get_vertices( self ):
        return arange( self.n_nodes_dof )[newaxis, :]
    def _get_vertex_numbers( self ):
        return arange( self.n_nodes_dof )[newaxis, :]

    def _get_lines( self ):
        # lines used for visualization
        print 'getting lines'
        lines_x = []
        lines_y = []
        lines_z = []
        enum_nodes_geo = self.enum_nodes_geo
        n_e_nodes_geo = self.n_e_nodes_geo
        # specify the limits of the loop in order to avoid an empty range
        range_x = range( max( 1, self.shape[0] ) )
        range_y = range( max( 1, self.shape[1] ) )
        range_z = range( max( 1, self.shape[2] ) )
        # loop over the elements
        for i in range_x:
            for j in range_y:
                for k in range_z:
                    # get the lines in x-, y-, and z-direction, respectively        
                    for m in range( 0, 2 ):
                        for n in range ( 0, 2 ):
                            # check if the system has an extension in the x-direction, i.e. lines in x-direction
                            if self.shape[0] != 0:
                                # check if the system has an extension in y- and z-direction, respectively
                                if self.shape[1] == 0:
                                    n = 0
                                if self.shape[2] == 0:
                                    m = 0
                                lines_x += [ [enum_nodes_geo[n_e_nodes_geo[0] * i     ,
                                                             n_e_nodes_geo[1] * ( j + n )  ,
                                                             n_e_nodes_geo[2] * ( k + m )  ],
                                              enum_nodes_geo[n_e_nodes_geo[0] * ( i + 1 )  ,
                                                             n_e_nodes_geo[1] * ( j + n )  ,
                                                             n_e_nodes_geo[2] * ( k + m )  ]]]
                            # check if the system has an extension in the y-direction, i.e. lines in y-direction
                            if self.shape[1] != 0:
                                # check if the system has an extension in x- and z-direction, respectively
                                if self.shape[0] == 0:
                                    n = 0
                                if self.shape[2] == 0:
                                    m = 0
                                lines_y += [ [enum_nodes_geo[n_e_nodes_geo[0] * ( i + n )  ,
                                                             n_e_nodes_geo[1] * j     ,
                                                             n_e_nodes_geo[2] * ( k + m )  ],
                                              enum_nodes_geo[n_e_nodes_geo[0] * ( i + n )  ,
                                                             n_e_nodes_geo[1] * ( j + 1 )  ,
                                                             n_e_nodes_geo[2] * ( k + m )  ]]]
                            # check if the system has an extension in the z-direction, i.e. lines in z-direction
                            if self.shape[2] != 0:
                                # check if the system has an extension in x- and y-direction, respectively
                                if self.shape[0] == 0:
                                    n = 0
                                if self.shape[1] == 0:
                                    m = 0
                                lines_z += [ [enum_nodes_geo[n_e_nodes_geo[0] * ( i + n )  ,
                                                             n_e_nodes_geo[1] * ( j + m )  ,
                                                             n_e_nodes_geo[2] * k     ],
                                              enum_nodes_geo[n_e_nodes_geo[0] * ( i + n )  ,
                                                             n_e_nodes_geo[1] * ( j + m )  ,
                                                             n_e_nodes_geo[2] * ( k + 1 )  ]]]
        lines_system = lines_x + lines_y + lines_z
        return lines_system

    def _get_faces( self ):
        faces = []
        for e in self.elements:
            n = len( e.nodes_geo )
            if n == 4:
                faces.append( list( e.nodes_geo ) )
        return faces

    def _get_volumes( self ):
        return array( [] )

    #---------------------------------------------------------------
    # View
    #---------------------------------------------------------------
    traits_view = View( Group( Item( 'lengths' ),
                               Item( 'shape' ),
                               Item( 'node_map_geo' ),
                               label = 'Discretization parameters' ),
#                       Group( Item('elements', style = 'custom', editor = elem_list_editor,
#                              show_label = False),
#                              label = 'Elements' ),
#                       Group( Item('nodes', style = 'custom', editor = node_list_editor,
#                                   show_label = False),
#                                   label = 'Nodes' ),
#                       Group( Item('X_mtx', style = 'custom', editor = ArrayViewEditor(),
#                                   show_label = False),
#                            label = 'Nodal coordinates' ),
                        Group( 
                              HSplit( 
                                     VSplit( Item( name = 'engine_view',
                                           style = 'custom',
                                           resizable = True,
                                           show_label = False,
                                           ),
                                           Item( name = 'current_selection',
                                           editor = InstanceEditor(),
                                           enabled_when = 'current_selection is not None',
                                           style = 'custom',
                                           springy = True,
                                           show_label = False ),
                                        ),
                                     Item( name = 'scene',
                                           editor = SceneEditor(),
                                           show_label = False,
                                           resizable = True,
                                           height = 500,
                                           width = 500 ),
                                           ),
                                    label = 'View' ),
                                    resizable = True,
                                    height = 0.8,
                                    width = 0.8,
                                    buttons = [OKButton, CancelButton],
                                    kind = 'subpanel',
                                    )


class RTraceDomainField( RTrace, MeshActors ):
    label = Str( 'RTraceDomainField' )
    var = Str( '' )
    var_eval = Callable
    idx = Int( -1 )
    sd = WeakRef( ISDomain )
    warp = Bool( False )
    warp_f = Float( 1. )

    # Set up the postprocessing geometry.  Currently, the mesh is taken as is
    # However, in a general case the postprocessing mesh is finer 
    # and depends on the order of shape functions within the element.
    # Therefore, elements specify their local triangulation in the form
    # field_entity_type   
    #

    def _get_points( self ):
        points = []
        dim_slice = self.var_eval.dim_slice
        for e in self.sd.elements:
            X = e.get_X_mtx()
            if dim_slice: X = X[:, dim_slice]
            points += list( self.var_eval.get_vtk_r_glb_arr( X ) )
        self.points = array( points )

        return self.points

    def _get_lines( self ):
        lines = []
        n_vtk_r = self.var_eval.n_vtk_r
        for i, e in enumerate( self.sd.elements ):
            lines += ( list( self.var_eval.field_lines + i * n_vtk_r ) )
        self.lines = array( lines )
        return self.lines

    def _get_faces( self ):
        faces = []
        n_vtk_r = self.var_eval.n_vtk_r
        for i, e in enumerate( self.sd.elements ):
            faces += ( list( self.var_eval.field_faces + i * n_vtk_r ) )
        self.faces = array( faces )
        return self.faces

    def bind( self ):
        '''
        Locate the evaluators
        '''
        self.var_eval = self.rmgr.rte_dict[self.var]

    def setup( self, sctx ):
        '''
        Setup the spatial domain of the tracer
        '''
        self.sctx = sctx
        self.sd = sctx.sdomain
        #
        # Find out which fets_evals are present in the domain and 
        # get the local coordinates of the sampled points. 
        #
        sd = sctx.sdomain
        self.elX_mtx = []
        for e in sd.elements:
            self.elX_mtx.append( e.get_X_mtx() )

    def add_current_values( self, sctx, U_k ):
        '''
        Invoke the evaluators in the current context for the specified control vector U_k.
        '''
        # Get the domain points
        # TODO - make this more compact. The element list is assumed to be uniform 
        # so that all element arrays have the same shape. Thus, use slices and vectorized 
        # evaluation to improve the performance 
        sd = self.sctx.sdomain
        loc_coords = self.var_eval.vtk_r_arr
        n_loc = loc_coords.shape[0]
        scalar_field = []
        for e in sd.elements:
            if self.var_eval.fets_eval.mats_eval:
                e_id = e.id_number
                e_arr_size = self.var_eval.get_state_array_size()
                self.sctx.elem_state_array = self.sctx.state_array[e_id * e_arr_size :\
                                                               ( e_id + 1 ) * e_arr_size]
            self.sctx.X = e.get_X_mtx()
            self.sctx.x = e.get_x_mtx()
            self.sctx.elem = e
            field_entry = []
            for i in range( n_loc ):
                if self.var_eval.fets_eval.mats_eval:
                    gp_id = self.var_eval.fets_eval.vtk_point_ip_map[i]
                    m_arr_size = self.var_eval.fets_eval.m_arr_size
                    self.sctx.mats_state_array = self.sctx.elem_state_array\
                                                [gp_id * m_arr_size: ( gp_id + 1 ) * m_arr_size]
                self.sctx.loc = loc_coords[i]
                val = self.var_eval( self.sctx, U_k )
                field_entry.append( val )
            scalar_field += field_entry
        self.field_arr = array( scalar_field )
        self.changed = True

    def _get_vectors( self ):
        return self.rmgr.warp_field

    field_entity_type = Delegate( 'var_eval' )
    def _get_scalars( self ):
        #return self.field_arr
        return self.field_arr[:, self.idx]

    def timer_tick( self, e = None ):
        pass

    def clear( self ):
        pass

#    view = View( HGroup ( VGroup( VGroup('var','idx'),
#                                  VGroup('record_on','clear_on') ),
#                                  Item('actor_dict',
#                                       editor=ActorEditor(scene_kwds={'background':(1.0,1.0,1.0)}),
#                                       show_label = False)),
#                                       resizable = True )

    view = View( HSplit( VSplit ( VGroup( 'var', 'idx' ),
                                  VGroup( 'record_on', 'clear_on' ),
                                  VSplit( Item( name = 'engine_view',
                                           style = 'custom',
                                           resizable = True,
                                           show_label = False
                                           ),
                                           Item( name = 'current_selection',
                                           editor = InstanceEditor(),
                                           enabled_when = 'current_selection is not None',
                                           style = 'custom',
                                           springy = True,
                                           show_label = False ),
                                           ) ),
                                    Item( name = 'scene',
                                           editor = SceneEditor(),
                                           show_label = False,
                                           resizable = True,
                                           height = 500,
                                           width = 500 ),
                                           ),
                                    resizable = True )

if __name__ == '__main__':
    # Discretization
    #

#    # show grid for default arguments:
#    grid_domain = MGridDomain()

    #grid_domain = MGridDomain( #field_entity_type = 'quad',
    #                              lengths = (3.,3.,3.),
    #                              shape = (3,2,1),
    #                              n_nodal_dofs = 2 )


    # generate grid based on specified arguments:
#    grid_domain = MGridDomain( lengths = (3.,3.,0.),
#                                  shape = (1,1,0),
#                                  n_nodal_dofs = 2 ,
#                                  n_e_nodes_dof = (1,1,0),
#                                  n_e_nodes_geo = (1,1,0),
#                                  node_map_geo = [0,1,3,2],
#                                  node_map_dof = [0,1,3,2] )

    # Discretization
    #
    grid_domain = MGridDomain( lengths = ( 3., 3., 0. ),
                             shape = ( 2, 1, 0 ),
                             n_nodal_dofs = 2,
                             # NOTE: the following properties must be defined and 
                             # must correspond to the used element formulation
                             n_e_nodes_geo = ( 1, 1, 0 ),
                             n_e_nodes_dof = ( 2, 2, 0 ),
                             node_map_geo = [0, 1, 3, 2],
                             node_map_dof = [0, 2, 8, 6, 1, 5, 7, 3] )

    #print grid_domain.elements
#
#    for i, e in enumerate( grid_domain.elements ):
#        print 'e.id_number', e.id_number
#        print 'e.nodes_geo', e.nodes_geo
#        print 'e.get_X_mtx()', e.get_X_mtx()


#    print 'top'
#    print grid_domain.get_top_dofs()
    #print grid_domain.elcoord_array
    print grid_domain.n_dofs

    grid_domain.changed = True
    grid_domain.configure_traits( view = "traits_view" )