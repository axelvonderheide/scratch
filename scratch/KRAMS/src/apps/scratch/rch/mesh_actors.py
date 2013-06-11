from enthought.traits.api import \
     Array, Bool, Enum, Float, HasTraits, HasStrictTraits, \
     Instance, Int, Trait, Str, \
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

from numpy import \
    array


from enthought.tvtk.api import tvtk
# Mayavi related imports
#
from enthought.tvtk.pyface.scene_editor import SceneEditor 
from enthought.mayavi.view.engine_view import EngineView
#from enthought.mayavi.core.ui.engine_view import EngineView
from enthought.mayavi import mlab
from enthought.mayavi.tools.mlab_scene_model import MlabSceneModel
from enthought.mayavi.sources.vtk_data_source import VTKDataSource
from enthought.mayavi.modules.api import Outline, Surface
from enthought.mayavi.filters.api import WarpVector

# Set mlab to use the simple backend instead of envisage.
mlab.options.backend = 'simple'

#-------------------------------------------------------------------
# FEMeshView - visualization of the domain
#-------------------------------------------------------------------

class MeshActors(HasTraits):
    '''
    Provide a tvtk actor with the mesh for visualization. 

    @TODO must be revised in the context of the view and editor
    concept.  It is now a mixin-base class providing the visualization
    functionality. Still its the question whether or not it should be
    implemented taht way.
    '''
    #---------------------------------------------------------------
    # Visualization methods
    #---------------------------------------------------------------
    # Convertion of mesh representation to tvtk-understandable arrays
    #
    def _get_points( self ):
        raise NotImplementedError
    def _get_point_numbers( self ):
        return array( [] )
    def _get_vertices( self ):
        return array( [] )
    def _get_vertex_numbers( self ):
        return array( [] )
    def _get_lines( self ):
        return array( [] )
    def _get_faces( self ):
        return array( [] )
    def _get_volumes( self ):
        return array( [] )
    def _get_field_arr( self ):
        return array( [] )

    field_entity_type = None

    scene = Instance(MlabSceneModel, ())
    
    # The mayavi engine view.
    engine_view = Instance(EngineView)

    # The current selection in the engine tree view.
    current_selection = Property
    
    def __init__(self, **traits):
        HasTraits.__init__(self, **traits)
        self.engine_view = EngineView(engine=self.scene.engine)

        # Hook up the current_selection to change when the one in the engine
        # changes.  This is probably unnecessary in Traits3 since you can show
        # the UI of a sub-object in T3.
        self.scene.engine.on_trait_change(self._selection_change,
                                          'current_selection')
        
    changed = Event     
    @on_trait_change('changed')
    def generate_scene( self ):
        '''
        Generates the Mayavi scene
        '''
        points,lines,faces = self.tvtk_source
        f_data = tvtk.PolyData(points = points, lines = lines, polys = faces)

        if self.field_entity_type == 'quad' or self.field_entity_type == 'line':
            arr = self._get_scalars()
            f_data.point_data.scalars = arr
            f_data.point_data.scalars.name = self.var+' '+str(self.idx)
            if self.label == 'RTraceDomainField' and self.warp:
                warp = self._get_vectors()
                f_data.point_data.vectors = warp
                f_data.point_data.vectors.name = 'displ'
            
        src = VTKDataSource(data = f_data)
        e = self.scene.engine
        e.add_source(src)
        if self.label == 'RTraceDomainField' and self.warp:
            e.add_filter(WarpVector())
        e.add_module(Surface())
        #e.add_module(Outline())
        
    tvtk_source = Property(Any)
    @cached_property
    def _get_tvtk_source( self ):
        # data
        points   = self._get_points()
        point_numbers   = self._get_point_numbers()
        vertices = self._get_vertices()
        lines    = self._get_lines()
        faces    = self._get_faces()
        volumes  = self._get_volumes()

        # Faces
        #
        
        #f_data = tvtk.UnstructuredGrid(points=points)
        #cel_type = 7
        #f_data.set_cells(cel_type, faces)
        return points, lines, faces
        
    
    def _selection_change(self, old, new):
        self.trait_property_changed('current_selection', old, new)
        
    
    def _get_current_selection(self):
        return self.scene.engine.current_selection