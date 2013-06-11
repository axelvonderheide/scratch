
from enthought.traits.api import \
     Array, Bool, Callable, Enum, Float, HasTraits, Interface, implements, \
     Instance, Int, Trait, Str, Enum, Callable, List, TraitDict, Any, \
     on_trait_change, Tuple, WeakRef, Delegate, Property, cached_property, Dict

from enthought.traits.ui.api import \
     Item, View

from enthought.traits.ui.menu import \
     OKButton, CancelButton
     

from numpy import \
     zeros, float_, ix_, array, frompyfunc, average, where,\
     index_exp, hstack, vstack, unique

from ibvpy.api import \
     ITStepperEval, TStepperEval

from ibvpy.core.rtrace_eval import RTraceEval
from ibvpy.fets.fets_eval import IFETSEval
from ibvpy.dots.dots_eval import DOTSEval
#-----------------------------------------------------------------------------
# Integrator for a simple 1D domain.
#-----------------------------------------------------------------------------

class DOTSEval_enr( DOTSEval ):
    '''
    Domain with uniform FE-time-step-eval.
    '''
    implements( ITStepperEval )

    fets_eval = Instance( IFETSEval )
    some = Str('')

    def new_cntl_var(self):
        return zeros( self.domain.n_dofs, float_ )
	
    def new_resp_var(self):
        return zeros( self.domain.n_dofs, float_ )
	
    def get_state_array_size( self ):

        # The overall size is just a n_elem times the size of a single element
        #
        shape = sctx.sdomain.shape
        #self.e_arr_size = shape * self.fets_eval.get_state_array_size( )
        self.e_arr_size = self.fets_eval.get_state_array_size( )
        dots_arr_size = shape * self.e_arr_size
        return dots_arr_size

    def setup( self, sctx ):
        self.domain = sctx.sdomain
        
        ndofs = self.domain.n_dofs
        self.K = zeros((ndofs, ndofs), float_ )
        self.F_int = zeros( ndofs, float_ )

        e_arr_size      = self.e_arr_size

        # Run the setup of sub-evaluator
        #
        self.fets_eval.gp_list = self.fets_eval.st_gp_list
        for e_id, elem in enumerate( sctx.sdomain.elements ):
            sctx.elem = elem
            sctx.elem_state_array = sctx.state_array[ e_id * e_arr_size : (e_id + 1) * e_arr_size ]
            self.fets_eval.setup( sctx )
	
    def get_corr_pred( self, sctx, u, du, tn, tn1 ):
	
        self.K[:,:] = 0.0
        self.F_int[:]    = 0.0
        e_arr_size      = self.e_arr_size
        intersect_elems = self.domain.get_intersected_elems()

        for e_id, elem in enumerate( sctx.sdomain.elements ):
            if e_id in intersect_elems:
                #TODO: effectivity - this is called in every iteration
                #TODO: deactivate "negative elems" -cutting with ls
                pt_sets = self._get_ls_points(sctx, self.fets_eval._node_coord_map,elem.get_X_mtx())
                gpl,vtk_r, field_faces =self.fets_eval.get_tri_gp_list(pt_sets)
                self.fets_eval.gp_list = gpl
                elem.vtk_r = vtk_r
                elem.field_faces = field_faces
            else:
                self.fets_eval.gp_list = self.fets_eval.st_gp_list
            ix = elem.get_dof_map()
            sctx.elem = elem
            sctx.elem_state_array = sctx.state_array[ e_id*e_arr_size : (e_id+1)*e_arr_size ]
            sctx.X = elem.get_X_mtx()
            f, k = self.fets_eval.get_corr_pred( sctx, u[ix_(ix)], du[ix_(ix)], tn, tn1 )
            self.K[ ix_(ix,ix) ] += k
            self.F_int[ ix_(ix) ] += f


        return self.F_int, self.K
    
    def _get_ls_points(self, sctx, r_pnt, X_mtx):
        #print "LOC GLOB ",r_pnt," ", X_mtx
        values = []
        for X_pnt in X_mtx:
            values.append(self.domain.get_ls_value(X_pnt))
        vals = array(values)
        pos_idx = where(vals > 0)
        neg_idx = where(vals < 0)
        values_ext = hstack((values,values[0]))
        r_pnt_ext =  vstack((r_pnt,r_pnt[0]))
        idx_arr = where(values_ext[:-1] * values_ext[1:] <= 0)
        coord = []
        for i in idx_arr[0]:
            int_coord = average([r_pnt_ext[i],r_pnt_ext[i+1]],0,\
                                [abs(values_ext[i+1]),abs(values_ext[i])])#weigths got to be swapped!
            #TODO: iteration to ensure point lies on the isoline 
#            sctx.r_pnt = int_coord
#            X_pnt = self.fets_eval.get_X_pnt(sctx)
#            error = self.domain.get_ls_value(X_pnt[0])
#            print "error before ",error
#            e_count = 0
#            while abs(error) > 1.e-2 and e_count < 10:
#                if (error * values_ext[i])<0:
#                    int_coord = average([r_pnt_ext[i],int_coord],0,\
#                                [abs(error),abs(values_ext[i])])
#                else:
#                    int_coord = average([r_pnt_ext[i+1],int_coord],0,\
#                                [abs(error),abs(values_ext[i+1])])
#                sctx.r_pnt = int_coord
#                X_pnt = self.fets_eval.get_X_pnt(sctx)
#                error = self.domain.get_ls_value(X_pnt[0])
#                e_count += 1
#            print "error after ", error
            coord.append(int_coord)
        pos_point = vstack((r_pnt[index_exp[pos_idx]],coord))
        neg_point = vstack((r_pnt[index_exp[neg_idx]],coord))
        #print "PN ",pos_point
        #print "NN ",neg_point
        return [pos_point, neg_point]
    

    def map_u(self, sctx, U):
        ix = sctx.elem.get_dof_map()
#        sctx.r = fets_eval.map_to_local( sctx.elem, sctx.X )
        u = U[ix]
        return u
    
    rte_dict = Property( Dict, depends_on = 'fets_eval')
    @cached_property
    def _get_rte_dict(self):
        rte_dict = {}
        for key, eval in self.fets_eval.rte_dict.items():
            rte_dict[ key ] = RTraceEvalUDomainFieldVar( name = key, 
                                                  u_mapping = self.map_u, 
                                                  eval = eval,
                                                  fets_eval = self.fets_eval )
        return rte_dict
      
    traits_view = View( Item('fets_eval', style = 'custom', show_label = False ),
                        resizable = True,
                        height = 0.8,
                        width = 0.8,
                        buttons = [OKButton,CancelButton],
                        kind = 'subpanel',
                        scrollable = True,
                        )

class RTraceEvalUDomainFieldVar(RTraceEval):
    fets_eval = WeakRef( IFETSEval )
    
    # @TODO Return the parametric coordinates of the element covering the element domain
    #
    vtk_r_arr = Delegate('fets_eval')
    n_vtk_r = Delegate('fets_eval')
    field_entity_type = Delegate('fets_eval')
    dim_slice = Delegate('fets_eval')
    get_vtk_r_glb_arr = Delegate('fets_eval')
    field_faces = Delegate('fets_eval')
    field_lines = Delegate('fets_eval')
    get_state_array_size = Delegate('fets_eval')
    
    # @TODO Return the entity lists of the element used for element visualization
    #
    def _get_points(self, X):
        return fets_eval.get_points()
    
from ibvpy.api import RTrace
from ibvpy.core.i_sdomain import \
    ISDomain
from ibvpy.plugins.mayavi.pipelines import \
    MVPolyData, MVPointLabels
    
from enthought.traits.ui.api import Item, View, HGroup, ListEditor, VGroup, \
     HSplit, Group, Handler, VSplit, TableEditor, ListEditor
    
class RTraceDomainField_enr(RTrace):
    '''
    Trace encompassing the whole spatial domain.
    '''
    label    = Str('RTraceDomainField')
    var      = Str('')
    var_eval = Callable
    idx      = Int(-1)
    sd       = WeakRef( ISDomain )
    warp     = Bool(False)
    warp_f   = Float(1.)


    #----------------------------------------------------------------------------
    # Visualization pipelines
    #----------------------------------------------------------------------------
    mvp_mgrid_geo = Trait( MVPolyData )
    def _mvp_mgrid_geo_default(self):
        return MVPolyData( name = 'Field %s' % self.var,
                               points = self._get_points,
                               #lines  = self._get_lines,
                               polys  = self._get_faces,
                               scalars = self._get_scalars,
                               #vectors = self._get_vectors
                               )
        
    mvp_imgrid_geo = Trait( MVPolyData )
    def _mvp_imgrid_geo_default(self):
        return MVPolyData( name = 'Field %s' % self.var,
                               points = self._get_ipoints,
                               #lines  = self._get_lines,
                               polys  = self._get_ifaces,
                               #scalars = self._get_scalars,
                               #vectors = self._get_vectors
                               )

    def redraw(self):
        '''
        '''
        self.mvp_mgrid_geo.redraw() # 'label_scalars')  
        self.mvp_imgrid_geo.redraw()      
        
    # Set up the postprocessing geometry.  Currently, the mesh is taken as is
    # However, in a general case the postprocessing mesh is finer 
    # and depends on the order of shape functions within the element.
    # Therefore, elements specify their local triangulation in the form
    # field_entity_type   
    #
    def _get_points(self):
        points = []
        ipoints = []
        dim_slice = self.var_eval.dim_slice
        intersect_elems = self.sctx.sdomain.get_intersected_elems()#elements intersected by levelset
        for i,e in enumerate(self.sd.elements):
            if i in intersect_elems:
                r_pnt = e.vtk_r
                X = e.get_X_mtx()
                if dim_slice: X = X[:,dim_slice]
                ipoints += list( self.var_eval.get_vtk_r_glb_arr( X,r_pnt ) )
            else:    
                X = e.get_X_mtx()
                if dim_slice: X = X[:,dim_slice]
                points += list( self.var_eval.get_vtk_r_glb_arr( X ) )
        self.points = array(points)
        self.ipoints = array(ipoints)
        return self.points
    
    def _get_ipoints(self):
        #self._get_points()
        #self._get_faces()
        return self.ipoints
    
    
    def _get_lines(self):
        lines = []
        n_vtk_r = self.var_eval.n_vtk_r
        for i, e in enumerate( self.sd.elements ):
            lines += ( list( self.var_eval.field_lines + i * n_vtk_r ) )
        self.lines = array(lines)         
        return self.lines
    
    def _get_faces(self):
        faces = []
        ifaces = []
        skipped = 0
        point_offset = 0
        n_vtk_r = self.var_eval.n_vtk_r
        intersect_elems = self.sctx.sdomain.get_intersected_elems()#elements intersected by levelset
        for i, e in enumerate( self.sd.elements ):
            if i in intersect_elems:
                ifaces += (list( e.field_faces +  point_offset) )
                point_offset += len(unique(e.field_faces))#Triangulation
                skipped += 1
            else:       
                #faces += ( list( self.var_eval.field_faces + i * n_vtk_r ) )
                faces_offset = (i-skipped) * n_vtk_r
                faces += ( list( self.var_eval.field_faces +  faces_offset) )
        self.faces = array(faces) 
        self.ifaces = array(ifaces)       
        return self.faces
    
    def _get_ifaces(self):
        print "faces ", self.ifaces
        return self.ifaces
        
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
        intersect_elems = sd.get_intersected_elems()#elements intersected by levelset
        
        scalar_field = []
        iscalar_field = []
        dim_slice = self.var_eval.dim_slice 
        for e_id, e in enumerate( sd.elements ):
            if e_id in intersect_elems:
                pass
            else:
                loc_coords = self.var_eval.vtk_r_arr
                n_loc = loc_coords.shape[0]
                e_arr_size = self.var_eval.get_state_array_size()
                self.sctx.elem_state_array = self.sctx.state_array[e_id*e_arr_size :\
                                                               (e_id+1)*e_arr_size]                
                X = e.get_X_mtx()
                if dim_slice:
                    X = X[:,dim_slice]
                self.sctx.X = X
                self.sctx.elem = e
                field_entry = []
                for i in range(n_loc):
                    gp_id = self.var_eval.fets_eval.vtk_point_ip_map[i]
                    m_arr_size = self.var_eval.fets_eval.m_arr_size
                    self.sctx.mats_state_array = self.sctx.elem_state_array\
                                                [gp_id * m_arr_size: (gp_id+1)*m_arr_size]
                    self.sctx.loc = loc_coords[i]
                    val = self.var_eval( self.sctx, U_k )
                    field_entry.append(  val )
                scalar_field += field_entry
        self.field_arr = array(scalar_field)

    def _get_vectors(self):
        return self.rmgr.warp_field

    field_entity_type = Delegate( 'var_eval' )
    def _get_scalars(self):
        return self.field_arr[:,self.idx]

    def timer_tick( self, e = None ):
        #self.changed = True
        pass

    def clear(self):
        pass

#    view = View( HGroup ( VGroup( VGroup('var','idx'),
#                                  VGroup('record_on','clear_on') ),
#                                  Item('actor_dict',
#                                       editor=ActorEditor(scene_kwds={'background':(1.0,1.0,1.0)}),
#                                       show_label = False)),
#                                       resizable = True )
        
    view = View( HSplit( VSplit ( VGroup('var','idx'),
                                  VGroup('record_on','clear_on'),
                                  VGroup( Item('refresh_button', show_label = False ) ),
                                                                    
#                                  VSplit(Item(name='engine_view',
#                                           style='custom',
#                                           resizable=True,
#                                           show_label=False
#                                           ),
#                                           Item(name='current_selection',
#                                           editor=InstanceEditor(),
#                                           enabled_when='current_selection is not None',
#                                           style='custom', 
#                                           springy=True,
#                                           show_label=False),
#                                           )),
#                                    Item(name='scene', 
#                                           editor=SceneEditor(),
#                                           show_label=False,
#                                           resizable=True,
#                                           height=500,
#                                           width=500
                                           ),
                                           ),                                  
                                    resizable = True )

