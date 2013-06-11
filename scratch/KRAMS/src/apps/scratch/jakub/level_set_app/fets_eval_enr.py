
from enthought.traits.api import \
     Array, Bool, Callable, Enum, Float, HasTraits, Interface, implements, \
     Instance, Int, Trait, Str, Enum, Callable, List, TraitDict, Any, \
     on_trait_change, Tuple, WeakRef, Delegate, Property, cached_property, Dict

from enthought.traits.ui.api import \
     View, Item
    
from numpy import \
     array, zeros, float_, dot, hstack, \
     average, ix_, int_, vstack, unique

from scipy.linalg import \
     det

from ibvpy.core.i_tstepper_eval import \
     ITStepperEval 
     
from ibvpy.core.tstepper_eval import \
     TStepperEval

from ibvpy.mats.mats_eval import \
     IMATSEval
     
from ibvpy.core.rtrace_eval import RTraceEval
from ibvpy.fets.fets_eval import FETSEval,IFETSEval
from enthought.tvtk.api import tvtk

class FETSEval_enr( FETSEval ):

    implements( IFETSEval )
    
    id_number = Int
    mats_eval = Instance( IMATSEval )
    # Distinguish the type of base geometric entity to be used for 
    # the visualization of the results.
    #
    field_entity_type = Enum('vertex','line','triangle','quad','tetra','hexa')
    vtk_r = Array( Float )
    field_points = Array( Int )
    field_lines = Array( Int )
    field_faces = Array( Int )
    field_volumens = Array( Int )
    
    n_vtk_r = Property(Int,depends_on = 'vtk_r')
    @cached_property
    def _get_n_vtk_r(self):
        return self.vtk_r.shape[0]
    
    gp_list = Trait()
    
    st_gp_list = Property(Array,depends_on = 'ngp_r,ngp_s,ngp_t')
    @cached_property
    def _get_st_gp_list(self):
        gp_list = []
        r_range = max(1,self.ngp_r)
        s_range = max(1,self.ngp_s)
        t_range = max(1,self.ngp_t)
        for gp_t in range(0,t_range):
            for gp_s in range(0,s_range):
                for gp_r in range(0,r_range): 
                    gp = [self.GP[self.ngp_r]['gp'][gp_r],\
                          self.GP[self.ngp_s]['gp'][gp_s],\
                          self.GP[self.ngp_t]['gp'][gp_t]],\
                          self.GP[self.ngp_r]['gw'][gp_r]*\
                          self.GP[self.ngp_s]['gw'][gp_s]*\
                          self.GP[self.ngp_t]['gw'][gp_t]       
                    gp_list.append(gp)
        return gp_list
    
    def get_tri_gp_list(self, pt_sets):
        #geo_r = self._node_coord_map
        points_list = []
        triangles_list = []
        point_offset = 0
        for pts in pt_sets:
            pt, tri = self._get_triangulation(pts)
            points_list += list(pt)
            triangles_list += list(tri+point_offset)
            point_offset += len(unique(tri))#Triangulation
        points = array(points_list)
        triangles = array(triangles_list)
        #print "GP: ",triangles
        int_order = 3
        gps=[]
        if int_order == 1:
            for id in triangles:
                gp=[average(points[ix_(id)],0),1.]
                #print "gp ",gp
                gps.append(gp)
        elif int_order == 2:    
            raise NotImplementedError
        elif int_order == 3:    
            weigths = array([[0.6,0.2,0.2],[0.2,0.6,0.2],[0.2,0.2,0.6]])
            for id in triangles:
                gps +=[average(points[ix_(id)],0),-0.5625],\
                    [average(points[ix_(id)],0, weigths[0]),0.52083333333333337],\
                    [average(points[ix_(id)],0, weigths[1]),0.52083333333333337],\
                    [average(points[ix_(id)],0, weigths[2]),0.52083333333333337]
    
                #gps.append(gp)
        elif int_order == 4:    
            raise NotImplementedError
        elif int_order == 5:    
            weigths = array([[0.0597158717, 0.4701420641, 0.4701420641],\
                             [0.4701420641, 0.0597158717, 0.4701420641],\
                             [0.4701420641, 0.4701420641, 0.0597158717],\
                             [0.7974269853, 0.1012865073, 0.1012865073],\
                             [0.1012865073, 0.7974269853, 0.1012865073],\
                             [0.1012865073, 0.1012865073, 0.7974269853]])
            for id in triangles:
                weigts_sum = False#for debug
                gps +=[average(points[ix_(id)],0),0.225],\
                    [average(points[ix_(id)],0, weigths[0],weigts_sum),0.1323941527],\
                    [average(points[ix_(id)],0, weigths[1],weigts_sum),0.1323941527],\
                    [average(points[ix_(id)],0, weigths[2],weigts_sum),0.1323941527],\
                    [average(points[ix_(id)],0, weigths[3],weigts_sum),0.1259391805],\
                    [average(points[ix_(id)],0, weigths[4],weigts_sum),0.1259391805],\
                    [average(points[ix_(id)],0, weigths[5],weigts_sum),0.1259391805]
    
                #gps.append(gp)
        else:    
            raise NotImplementedError
        #print "gps ",gps
        return gps, points[:,:-1], triangles #points for 2d
    
    def _get_triangulation(self, points):
        n_dims = points.shape[1]
        n_add = 3 - n_dims
        if n_add > 0:
            points = hstack( [points, 
                        zeros( [points.shape[0],n_add], dtype = 'float_' )] )
        # Create a polydata with the points we just created.
        profile = tvtk.PolyData(points=points)


        # Perform a 2D Delaunay triangulation on them.
        delny = tvtk.Delaunay2D(input=profile)
        tri= delny.output
        tri.update()
        triangles = array(tri.polys.data,dtype=int_)    
        #pts = array(tri.points.data)
        ids = (triangles.reshape((triangles.shape[0]/4),4))[:,1:]
        return tri.points.data, ids 
    
    
    # The user-specified fv_loc_coords list gets transform to an internal 
    # array representation
    #
    vtk_r_arr = Property( depends_on = 'vtk_r' )
    @cached_property
    def _get_vtk_r_arr(self):
        return array( self.vtk_r )
    
    def get_vtk_r_glb_arr(self, X_mtx, r_mtx = None ):
        '''
        Get an array with global coordinates of the element decomposition.
        
        If the local_point_list is non-empty then use it instead of the one supplied 
        by the element specification. This is useful for augmented specification of RTraceEval 
        evaluations with a specific profile of a field variable to be traced.
        '''
        if self.dim_slice:
            X_mtx = X_mtx[:,self.dim_slice]

        if r_mtx == None: r_mtx = self.vtk_r_arr
        
        # TODO - efficiency in the extraction of the global coordinates. Is broadcasting 
        # a possibility - we only need to augment the matrix with zero coordinates in the
        # in the unhandled dimensions. 
        #
        X3D = array( [ dot( self.get_N_geo_mtx(r_pnt), X_mtx )[0,:] for r_pnt in r_mtx ] )
        n_dims = r_mtx.shape[1]
        n_add = 3 - n_dims
        if n_add > 0:
            X3D = hstack( [X3D, 
                           zeros( [r_mtx.shape[0],n_add], dtype = 'float_' )] )
        return X3D
    
    def get_X_pnt(self, sctx):
        '''
        Get the global coordinates for the specified local coordinats r_pnt
        @param r_pnt: local coordinates
        '''
        r_pnt = sctx.r_pnt
        X_mtx = sctx.X
        if self.dim_slice:
            X_mtx = X_mtx[:,self.dim_slice]

        # TODO - efficiency in the extraction of the global coordinates. Is broadcasting 
        # a possibility - we only need to augment the matrix with zero coordinates in the
        # in the unhandled dimensions. 
        #
        return dot( self.get_N_geo_mtx(r_pnt), X_mtx )
        
    # Number of element DOFs
    #
    n_e_dofs = Int
    
    # Dimensionality
    dim_slice = None
    
    # Parameters for the time-loop
    #
    def new_cntl_var(self):
        return zeros( self.n_e_dofs, float_ )
    
    def new_resp_var(self):
        return zeros( self.n_e_dofs, float_ )

    def get_state_array_size( self ):
        r_range = max(1,self.ngp_r)
        s_range = max(1,self.ngp_s)
        t_range = max(1,self.ngp_t)
        self.m_arr_size = self.mats_eval.get_state_array_size()
        return self.m_arr_size * r_range * s_range * t_range
        

    def setup(self, sctx):
        i = 0
        for gp in self.gp_list:
            sctx.mats_state_array = sctx.elem_state_array[(i * self.m_arr_size): ((i+1)*self.m_arr_size)]
            self.mats_eval.setup( sctx )
            i += 1

    ngp_r = Int(0,label = 'Number of Gauss points in r-direction')
    ngp_s = Int(0,label = 'Number of Gauss points in s-direction')
    ngp_t = Int(0,label = 'Number of Gauss points in t-direction')
    #-------------------------------------------------------------------
    # Overloadable methods
    #-------------------------------------------------------------------
    def get_corr_pred(self, sctx, u, du, tn, tn1, u_avg = None ):
        '''
        Corrector and predictor evaluation.

        @param u current element displacement vector
        '''
        if self.dim_slice:
            X_mtx = sctx.X[:,self.dim_slice]
        else:            
            X_mtx = sctx.X
            
            

        ### Use for Jacobi Transformation
        
        n_e_dofs = self.n_e_dofs
        K = zeros((n_e_dofs,n_e_dofs))
        F = zeros(n_e_dofs)
        sctx.fets_eval = self  
#        i = 0      
#        for gp_t in range(0,ngp_t):
#            for gp_s in range(0,ngp_s):
#                for gp_r in range(0,ngp_r):
#                    r_pnt = array([gp_c_r[0,gp_r],gp_c_s[0,gp_s],gp_c_t[0,gp_t]])
#                    J_det = self._get_J_det( r_pnt, X_mtx )
#                    B_mtx = self.get_B_mtx( r_pnt, X_mtx )
#                    eps_mtx = dot( B_mtx, u )
#                    if self.mats_eval:
#                        sctx.mats_state_array = sctx.elem_state_array[i * self.m_arr_size: (i+1)*self.m_arr_size]
#                    sig_mtx, D_mtx = self.get_mtrl_corr_pred( sctx, eps_mtx )
#                    K += dot( B_mtx.T, dot( D_mtx,  B_mtx ) )\
#                     * gp_w_r[0,gp_r]* gp_w_s[0,gp_s]* gp_w_t[0,gp_t]
#                    F += dot( B_mtx.T, sig_mtx )\
#                     * gp_w_r[0,gp_r]* gp_w_s[0,gp_s]* gp_w_t[0,gp_t]
#                    i += 1
        i = 0      
        for gp in self.gp_list:
            r_pnt = gp[0]
            sctx.r_pnt = r_pnt
            J_det = self._get_J_det( r_pnt, X_mtx )
            B_mtx = self.get_B_mtx( r_pnt, X_mtx )
            eps_mtx = dot( B_mtx, u )
            d_eps_mtx = dot( B_mtx, du )
            sctx.mats_state_array = sctx.elem_state_array[i * self.m_arr_size: (i+1)*self.m_arr_size]
            if u_avg != None:
                eps_avg = dot( B_mtx, u_avg )
                sig_mtx, D_mtx = self.get_mtrl_corr_pred( sctx, eps_mtx, d_eps_mtx, tn, tn1, eps_avg )
            else:
                sig_mtx, D_mtx = self.get_mtrl_corr_pred( sctx, eps_mtx, d_eps_mtx, tn, tn1 )
            K += dot( B_mtx.T, dot( D_mtx,  B_mtx ) ) * gp[1]
            F += dot( B_mtx.T, sig_mtx ) * gp[1]
            i += 1
        
        F = F * J_det
        K = K * J_det
        return F, K

    #-------------------------------------------------------------------
    # Standard evaluation methods
    #-------------------------------------------------------------------
    def get_J_mtx(self,r_pnt,X_mtx):
        return dot( self.get_dNr_geo_mtx(r_pnt), X_mtx )

    #-------------------------------------------------------------------
    # Required methods
    #-------------------------------------------------------------------

    def get_N_geo_mtx(r_pnt):
        raise NotImplementedError

    def get_dNr_geo_mtx(r_pnt):
        raise NotImplementedError

    def get_N_mtx(r_pnt):
        raise NotImplementedError

    def get_B_mtx( self, r_pnt, X_mtx ):
        '''
        Get the matrix for kinematic mapping between displacements and strains.
        @param r local position within the element.
        @param X nodal coordinates of the element.

        @TODO[jakub] generalize
        '''
        raise NotImplementedError

    def get_mtrl_corr_pred(self, sctx, eps_mtx, d_eps, tn, tn1, eps_avg = None):
        if eps_avg != None:
            sig_mtx, D_mtx = self.mats_eval.get_corr_pred(sctx, eps_mtx, d_eps, tn, tn1, eps_avg)
        else:
            sig_mtx, D_mtx = self.mats_eval.get_corr_pred(sctx, eps_mtx, d_eps, tn, tn1,)
        return sig_mtx, D_mtx
            
    #-------------------------------------------------------------------
    # Private methods
    #-------------------------------------------------------------------

    def _get_J_det(self,r_pnt3d,X_mtx):
        if self.dim_slice:
            r_pnt = r_pnt3d[self.dim_slice]
        return det( self.get_J_mtx(r_pnt,X_mtx) )

    GP = { # if no gauss point is defined in one direction (e.g. for ngp_t=0 for a 2D-problem)
           # then the default value for ngp_t=0 is used and a weighting coefficient of value 1.
           # In 'get_gp' the maximum of npg and 1 is used as range in the loop which leads to 
           # a simple multiplication with 1 for that direction. The coordinate of the gauss point 
           # in this direction is set to zero ('gp' : [0.]).
           0 : {'gp' : [0.],
                'gw' : [1.]},
           #     
           1 : {'gp' : [0.],
                'gw' : [2.]},
           2 : {'gp' : [-0.57735026918962584,0.57735026918962584],
                'gw' : [1.,1.]},
           3 : {'gp' : [-0.7745966692414834,0.,0.7745966692414834],
                'gw' : [0.55555555555555558,0.88888888888888884,0.55555555555555558]}
           }
    
#    GP = { 1 : {'gp' : array([[0.]]),
#                'gw' : array([[2.]])},
#           2 : {'gp' : array([[-0.57735026918962584,0.57735026918962584]]),
#                'gw' : array([[1.,1.]])},
#           3 : {'gp' : array([[-0.7745966692414834,0.,0.7745966692414834]]),
#                'gw' : array([[0.55555555555555558,0.88888888888888884,0.55555555555555558]])}
#           }

    def get_eps(self, sctx, u):
        X_mtx = sctx.X
        r_pnt = sctx.loc
        B_mtx = self.get_B_mtx(r_pnt, X_mtx)
        eps = dot( B_mtx, u )
        return eps 

    def get_u(self, sctx, u):
        N_mtx = self.get_N_mtx( sctx.loc )
        return dot( N_mtx, u )

    debug_on = Bool(False)
    def _debug_rte_dict(self):
        '''
        RTraceEval dictionary with field variables used to verify the element implementation
        '''        
        if self.debug_on:
            return {'Ngeo_mtx' : RTraceEvalElemFieldVar( eval = lambda sctx, u: self.get_N_geo_mtx( sctx.loc ),
                                                ts = self ),
                    'N_mtx': RTraceEvalElemFieldVar( eval = lambda sctx, u: self.get_N_mtx( sctx.loc )[0],
                                                ts = self ),
                    'B_mtx0': RTraceEvalElemFieldVar( eval = lambda sctx, u: self.get_B_mtx( sctx.loc, sctx.X )[0],
                                                ts = self ),
                    'B_mtx1': RTraceEvalElemFieldVar( eval = lambda sctx, u: self.get_B_mtx( sctx.loc, sctx.X )[1],
                                                ts = self ),
                    'B_mtx2': RTraceEvalElemFieldVar( eval = lambda sctx, u: self.get_B_mtx( sctx.loc, sctx.X )[2],
                                                ts = self ),
                    'J_det'  : RTraceEvalElemFieldVar( eval = lambda sctx, u: 
                                                array( [det( self.get_J_mtx( sctx.loc, sctx.X ) ) ] ),
                                                              ts = self ) }
        else:
            return {}

    # List of mats that are to be chained
    #
 
    
    # Declare and fill-in the rte_dict - it is used by the clients to
    # assemble all the available time-steppers.
    #
    rte_dict = Trait( Dict )
    def _rte_dict_default(self):
        '''
        RTraceEval dictionary with standard field variables.
        '''
        rte_dict = self._debug_rte_dict()
        rte_dict.update( {'eps' : RTraceEvalElemFieldVar( eval = self.get_eps, ts = self ), 
                          'u'   : RTraceEvalElemFieldVar( eval = self.get_u,   ts = self )} )

        for key, v_eval in self.mats_eval.rte_dict.items():
            rte_dict[ key ] = RTraceEvalElemFieldVar( name = key, 
                                               u_mapping = self.get_eps, 
                                               eval = v_eval )
        return rte_dict
    
    traits_view = View( Item("_node_coord_map@"),
                         Item('field_entity_type'),
                         Item('field_faces'),
                         Item('field_lines'),
                         Item('vtk_r'),
                         Item('mats_eval'),
                         Item('n_e_dofs'),
                         Item('n_vtk_r'),
                         Item('n_nodal_dofs'),
#                         Item('rte_dict'),
                        resizable = True,
                        scrollable = True
                        )

class RTraceEvalElemFieldVar( RTraceEval ):
    
    # To be specialized for element level
    #
    field_entity_type = Delegate('ts')
    vtk_r_arr = Delegate('ts')
    get_vtk_r_glb_arr = Delegate('ts')
    field_vertexes = Delegate('ts')
    field_lines    = Delegate('ts')
    field_faces    = Delegate('ts')
    field_volumes  = Delegate('ts')
