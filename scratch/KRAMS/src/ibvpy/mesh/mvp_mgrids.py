from enthought.traits.api import \
     Array, Bool, Enum, Float, HasTraits, HasStrictTraits, \
     Instance, Int, Trait, Str, \
     Callable, List, TraitDict, Any, Range, \
     Delegate, Event, on_trait_change, Button, \
     Interface, WeakRef, implements, Property, cached_property, Tuple, \
     Dict
from enthought.traits.ui.api import Item, View, HGroup, ListEditor, VGroup, \
     HSplit, Group, Handler, VSplit, TableEditor, ListEditor

from numpy import \
    array

from enthought.tvtk.api import tvtk
# Mayavi related imports
#
from enthought.mayavi.sources.vtk_data_source import VTKDataSource
from enthought.mayavi.modules.api import Outline, Surface
from enthought.mayavi.core.source import Source
from enthought.mayavi.filters.api import WarpVector

# MayaVI engine used for the pipeline construction 
#
from ibvpy.plugins.mayavi_engine import get_engine

class MVPStructuredGrid( HasTraits ):
    dims = Callable
    def _dims(self):
        return lambda : (1,1,1)

    points  = Callable
    def _points_default(self):
        return lambda : array([])
    
    pd = Instance( tvtk.StructuredGrid )
    def _pd_default(self):
        return tvtk.StructuredGrid()
    
    name = Str('')

    def __init__(self,**kw):

        super( MVPStructuredGrid, self ).__init__(**kw)
        e = get_engine()

        from enthought.mayavi.modules.api import \
        Outline, Surface, StructuredGridOutline, GridPlane

        self.src = VTKDataSource( name = self.name, data = self.pd )
        e.add_source(self.src)
        
        o = StructuredGridOutline()
        e.add_module(o)
        
        for axis in ['x','y','z']:
            g = GridPlane( name = '%s - grid plane' % axis )
            g.grid_plane.axis = axis
            e.add_module(g)
                      
            
    def redraw(self):

        self.pd.dimensions = self.dims()
        self.pd.points = self.points()  
        self.src.data_changed = True


class MVPBase( HasTraits ):
    '''
    Mayavi Pipeline Base class
    '''

    points  = Callable
    def _points_default(self):
        return lambda : array([])
    
    lines   = Callable
    def _lines_default(self):
        return lambda : array([])

    polys   = Callable
    def _polys_default(self):
        return lambda : array([])

    scalars = Callable

    vectors = Callable

    tensors = Callable
    
    pd = Instance( tvtk.PolyData )
    def _pd_default(self):
        return tvtk.PolyData()
    
    def redraw(self):

        self.pd.points = self.points()
        self.pd.lines = self.lines()
        self.pd.polys = self.polys()
        if self.scalars:
            self.pd.point_data.scalars = self.scalars()
        if self.vectors:
            self.pd.point_data.vectors = self.vectors()
        if self.tensors:
            self.pd.point_data.tensors = self.tensors()
        self.src.data_changed = True

#-------------------------------------------------------------------
# MVMeshSource - mayavi mesh source
#-------------------------------------------------------------------

class MVPMeshGridGeo( MVPBase ):
    '''
    Provide a Mayavi source for polar visualization visualization. 

    @TODO must be revised in the context of the view and editor
    concept.  It is now a mixin-base class providing the visualization
    functionality. Still its the question whether or not it should be
    implemented taht way.
    '''

    name = Str('')

    def __init__(self,**kw):

        e = get_engine()

        super( MVPMeshGridGeo, self ).__init__(**kw)
        from enthought.mayavi.modules.api import Outline, Surface, Labels

        self.src = VTKDataSource( name = self.name, data = self.pd )
        e.add_source(self.src)
        
        o = Outline()
        e.add_module(o)
        s = Surface()
        e.add_module(s)
        
class MVPMeshGridLabels( MVPBase ):

    def __init__(self,**kw):

        e = get_engine()

        super( MVPMeshGridLabels, self ).__init__(**kw)
        from enthought.mayavi.modules.api import Outline, Surface, Labels

        self.src = VTKDataSource( name = self.name, data = self.pd )
        e.add_source(self.src)
        
        self.labels = Labels( name = 'Node numbers', object = self.src, 
                         label_format = '%g',
                         number_of_labels = 100 )
        e.add_module(self.labels)

    def redraw(self, label_mode = 'label_ids' ):
        super(MVPMeshGridLabels,self).redraw()
        self.labels.mapper.label_mode = label_mode
