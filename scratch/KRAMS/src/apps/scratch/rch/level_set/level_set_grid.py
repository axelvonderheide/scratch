
    
from enthought.traits.api import \
    HasTraits, List, Array, Property, cached_property, \
    Instance, Trait, Button, on_trait_change, Tuple, \
    Int, Float, implements, Delegate

from enthought.traits.ui.api import \
    View, Item

from ibvpy.core.i_sdomain import \
    ISDomain
    
from ibvpy.core.sdomain import \
    SDomain

from numpy import \
    array, unique, min, max, mgrid, ogrid, c_, alltrue, repeat, ix_, \
    arange, ones, zeros, multiply, sort, index_exp, indices, add, hstack, \
    frompyfunc, where

from ibvpy.plugins.mayavi.pipelines import \
    MVPolyData, MVPointLabels, MVStructuredGrid

from ibvpy.mesh.cell_grid.cell_spec import CellSpec, GridCell
from ibvpy.mesh.cell_grid.cell_array import CellView, ICellView, CellArray, ICellArraySource

from math import sin

from ibvpy.mesh.cell_grid.cell_grid import CellGrid


class LevelSetGrid( CellGrid ):
    '''Evaluate a level set function
    '''
    
    a = Float(1.5, enter_set = True, auto_set = False )
    b = Float(2.0, enter_set = True, auto_set = False )
    shape = [300,100]
    coord_max = [10,3,0]
    
    grid_cell_spec = CellSpec( node_coords = [[-1,-1],[1,-1],[1,1],[-1,1]] )
        
    def level_set_fn(self, x, y):
        '''Level set function evaluation.
        '''
        return y - ( sin( self.b * x ) + self.a )
        return x - self.coord_max[0] / 2. - y / 2.
        return y - self.coord_max[1] / 2.

    level_set_grid = Property( Array, depends_on = 'a,b,shape,coord_max' )
    def _get_level_set_grid(self):
        # evaluate the level set function on the element point grid
        # @todo - change - the point grid contains also internal nodes 
        # that should not be involved in the evaluation. 
        X, Y = self.point_grid
        vect_fn = frompyfunc( self.level_set_fn, 2, 1 )
        values = vect_fn( X, Y )
        return array( values, dtype = 'float_' )
        
    def _get_transiting_edges(self):
        
        ls = self.level_set_grid
        x_edges = where( ls[:-1,:] * ls[1:,:] <= 0 )
        y_edges = where( ls[:,:-1] * ls[:,1:] <= 0 )

        ii,jj = x_edges
        # Get element numbers for each dimension separately
        # for each entry in x_edges 
        e_idx = []
        for i,j in zip(ii,jj):
            if j < self.shape[1]:
                e_idx.append( [i,j] )
            if j > 0:
                e_idx.append( [i,j-1] )

        ii,jj = y_edges
        for i,j in zip(ii,jj):
            if i < self.shape[0]:
                e_idx.append( [i,j] )
            if i > 0:
                e_idx.append( [i-1,j] )

        e_exp = array( e_idx, dtype = int ).transpose()
        return ( e_exp[0,:], e_exp[1,:] )
        
    def get_elem_intersection(self):
        
        e_idx = mgd._get_transiting_edges()
        return unique( mgd.get_elem_grid()[ e_idx ] )
        
    def get_elem_grid(self):
        shape = reduce( lambda x,y: x*y, self.shape )
        return arange( shape ).reshape( self.shape )
        
    def get_mvscalars(self):
        return self.level_set_grid.swapaxes(0,self.n_dims-1).flatten()

    def _get_ielem_points(self):
        icells = self.get_elem_intersection()
        mvpoints = []
        for cell_idx in icells:
            mvpoints += list( self.get_cell_mvpoints(cell_idx) ) 
        return array( mvpoints, dtype = 'float_' )
        
    def _get_ielem_polys(self):
        ncells = len( self.get_elem_intersection() )
        return arange( ncells * 4 ).reshape( ncells, 4 )
        
    #-----------------------------------------------------------------
    # Visualization-related methods
    #-----------------------------------------------------------------

    mvp_point_grid = Trait( MVStructuredGrid )
    def _mvp_point_grid_default(self):
        return MVStructuredGrid( name = 'Point grid', 
                                  dims = self._get_mvpoints_grid_shape,
                                  points = self._get_mvpoints,
                                  scalars = self.get_mvscalars )

    mvp_intersect_elems = Trait( MVPolyData )
    def _mvp_intersect_elems_default(self):
        return MVPolyData( name = 'Intersected elements', 
                                  points = self._get_ielem_points,
                                  polys = self._get_ielem_polys )
    
    refresh_button = Button('Draw')
    @on_trait_change('refresh_button')
    def redraw(self):
        '''Redraw the point grid.
        '''
        self.mvp_point_grid.redraw( )        
        self.mvp_intersect_elems.redraw()

    #------------------------------------------------------------------
    # UI - related methods
    #------------------------------------------------------------------
    traits_view = View(Item('grid_cell_spec'),
                       Item('shape@'),
                       Item('coord_min'),
                       Item('coord_max'),
                       Item('refresh_button'),
                       Item('cell_array'),
                       resizable = True,
                       scrollable = True,
                       height = 0.5,
                       width = 0.5)        

    
if __name__ == '__main__':
    mgd = LevelSetGrid()
    
    from ibvpy.plugins.ibvpy_app import IBVPyApp
    ibvpy_app = IBVPyApp( ibv_resource = mgd )
    ibvpy_app.main()
    
            