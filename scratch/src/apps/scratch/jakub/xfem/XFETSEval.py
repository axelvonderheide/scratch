from enthought.traits.api import \
     Array, Bool, Callable, Enum, Float, HasTraits, Interface, implements, \
     Instance, Int, Trait, Str, Enum, Callable, List, TraitDict, Any, \
     on_trait_change, Tuple, WeakRef, Delegate, Property, cached_property
     
from numpy import \
     array, zeros, int_, float_, ix_, dot, linspace, hstack, vstack, arange, \
     identity, unique, average, frompyfunc, abs, linalg

from scipy.linalg import \
     inv

from enthought.tvtk.api import tvtk

from ibvpy.fets.fets_eval import FETSEval

class XFETSEval( FETSEval ):
    
    
    parent_elem = Instance(FETSEval)
    
    parent_nnodes = Property(Int,depends_on = 'parent_elem')
    @cached_property
    def _get_parent_nnodes(self):
        return len(self.parent_elem.dof_r)
    
    n_e_dofs = Property(Int,depends_on = 'parent_elem')
    @cached_property
    def _get_n_e_dofs(self):
        return self.parent_elem.n_e_dofs * 2
    
    
    dof_r = Property(List,depends_on = 'parent_elem')
    @cached_property
    def _get_dof_r(self):
        return self.parent_elem.dof_r
    
    geo_r = Property(List,depends_on = 'parent_elem')
    @cached_property
    def _get_geo_r(self):
        return self.parent_elem.geo_r
    
    vtk_r = Property(depends_on = 'triangles')
    @cached_property
    def _get_vtk_r(self):
        coords = self.pp_triangles[0]
        if self.dim_slice: coords = coords[:,self.dim_slice]
        return coords
    
    
    vtk_cells = Property(depends_on = 'triangles')
    @cached_property
    def _get_vtk_cells(self):
        return self.pp_triangles[1].tolist()
    
    
    int_order = Int(5)
    
    int_triangles = Property(depends_on = 'pt_sets')
    @cached_property
    def _get_int_triangles(self):
        point_set = [ vstack((self.pt_sets[0],self.pt_sets[2])),
                     vstack((self.pt_sets[1],self.pt_sets[2]))]
        return self._get_triangulation(point_set)
        
    pp_triangles = Property(depends_on = 'pt_sets')
    @cached_property
    def _get_pp_triangles(self):
        dir_vct = self.pt_sets[2][1] - self.pt_sets[2][0]
        norm_vct = array([-dir_vct[1],dir_vct[0]])
        shift = norm_vct * 0.0001 / linalg.norm(norm_vct)
        #print "shift ",shift
        #print "dir_vct ", dir_vct
        #print "norm_vct ", norm_vct
        if dot(self.pt_sets[0][0],norm_vct)> 0.:
            pos_pts =  self.pt_sets[2] + shift
            neg_pts =  self.pt_sets[2] - shift
        else:
            pos_pts =  self.pt_sets[2] - shift
            neg_pts =  self.pt_sets[2] + shift
        
        point_set = [ vstack((self.pt_sets[0],pos_pts)),
                     vstack((self.pt_sets[1],neg_pts))]
        return self._get_triangulation(point_set)
        
        
    def _get_triangulation(self, point_set):
        #n_dims = points.shape[1]
        #n_add = 3 - n_dims
        n_add = 1#hack for 2D
        points_list = []
        triangles_list = []
        point_offset = 0
        for pts in point_set:
            if n_add > 0:
                points = hstack( [pts, 
                                  zeros( [pts.shape[0],n_add], dtype = 'float_' )] )
            # Create a polydata with the points we just created.
            profile = tvtk.PolyData(points=points)
            
            # Perform a 2D Delaunay triangulation on them.
            delny = tvtk.Delaunay2D(input=profile)
            tri= delny.output
            tri.update()#initiate triangulation
            triangles = array(tri.polys.data,dtype=int_)    
            pt = tri.points.data
            tri = (triangles.reshape((triangles.shape[0]/4),4))[:,1:]
            points_list += list(pt)
            triangles_list += list(tri+point_offset)
            point_offset += len(unique(tri))#Triangulation
        points = array(points_list)
        triangles = array(triangles_list)
        return [points, triangles]
        
    ip_coords = Property( depends_on = 'int_order')
    @cached_property
    def _get_ip_coords(self):
        '''Get the array of integration points'''
        points,triangles = self.int_triangles
        gps=[]
        if self.int_order == 1:
            for id in triangles:
                gp=average(points[ix_(id)],0)
                #print "gp ",gp
                gps.append(gp)
        elif self.int_order == 2:    
            raise NotImplementedError
        elif self.int_order == 3:    
            weigths = array([[0.6,0.2,0.2],[0.2,0.6,0.2],[0.2,0.2,0.6]])
            for id in triangles:
                gps +=average(points[ix_(id)],0),\
                    average(points[ix_(id)],0, weigths[0]),\
                    average(points[ix_(id)],0, weigths[1]),\
                    average(points[ix_(id)],0, weigths[2])
    
        elif self.int_order == 4:    
            raise NotImplementedError
        elif self.int_order == 5:    
            weigths = array([[0.0597158717, 0.4701420641, 0.4701420641],\
                             [0.4701420641, 0.0597158717, 0.4701420641],\
                             [0.4701420641, 0.4701420641, 0.0597158717],\
                             [0.7974269853, 0.1012865073, 0.1012865073],\
                             [0.1012865073, 0.7974269853, 0.1012865073],\
                             [0.1012865073, 0.1012865073, 0.7974269853]])
            for id in triangles:
                weigts_sum = False#for debug
                gps += average(points[ix_(id)],0),\
                     average(points[ix_(id)],0, weigths[0],weigts_sum),\
                     average(points[ix_(id)],0, weigths[1],weigts_sum),\
                     average(points[ix_(id)],0, weigths[2],weigts_sum),\
                     average(points[ix_(id)],0, weigths[3],weigts_sum),\
                     average(points[ix_(id)],0, weigths[4],weigts_sum),\
                     average(points[ix_(id)],0, weigths[5],weigts_sum)
        else:    
            raise NotImplementedError
        #print "gps ",gps
        #return gps, points[:,:-1], triangles #points for 2d
        return array( gps, dtype = 'float_' )
            
    ip_weights = Property( depends_on = 'int_order')
    @cached_property
    def _get_ip_weights(self):
        '''Get the array of integration points'''
        points,triangles = self.int_triangles
        gps=[]
        if self.int_order == 1:
            for id in triangles:
                gp= 1.
                #print "gp ",gp
                gps.append(gp)
        elif self.int_order == 2:    
            raise NotImplementedError
        elif self.int_order == 3:    
            for id in triangles:
                gps += -0.5625,0.52083333333333337,0.52083333333333337,\
                       0.52083333333333337
    
        elif self.int_order == 4:    
            raise NotImplementedError
        elif self.int_order == 5:    
            for id in triangles:
                weigts_sum = False#for debug
                gps += 0.225,0.1323941527,0.1323941527,0.1323941527,\
                       0.1259391805,0.1259391805,0.1259391805
        else:    
            raise NotImplementedError
        #print "gps ",gps
        #return gps, points[:,:-1], triangles #points for 2d
        return array( gps, dtype = 'float_' )
    
    def get_B_mtx( self, r_pnt, X_mtx ):
        #print "len ",len(self.parent_elem.dof_r)
        J_mtx = self.get_J_mtx(r_pnt,X_mtx)
        #print "loc_coords ", r_pnt
        dNr_mtx = self.get_parent_dNr_mtx( r_pnt )
        dNr_psi_mtx = self.get_dNr_psi_mtx( r_pnt )
        dNx_mtx = dot( inv( J_mtx ), dNr_mtx  )
        dNx_psi_mtx = dot( inv( J_mtx ), dNr_psi_mtx  )
        Bx_mtx = zeros( (3, self.parent_nnodes *2*2 ), dtype = 'float_' )#TODO: just 2D
        for i in range(0,self.parent_nnodes):#TODO: just 2D
            Bx_mtx[0,i*4]   = dNx_mtx[0,i]
            Bx_mtx[1,i*4+1] = dNx_mtx[1,i]
            Bx_mtx[2,i*4]   = dNx_mtx[1,i]
            Bx_mtx[2,i*4+1] = dNx_mtx[0,i]
            
            Bx_mtx[0,i*4+2] = dNx_psi_mtx[0,i]
            Bx_mtx[1,i*4+3] = dNx_psi_mtx[1,i]
            Bx_mtx[2,i*4+2] = dNx_psi_mtx[1,i]
            Bx_mtx[2,i*4+3] = dNx_psi_mtx[0,i]           
        return Bx_mtx
    
    def get_parent_N_mtx(self,r_pnt):
        return self.parent_elem.get_N_mtx(r_pnt)
    
    def get_parent_dNr_mtx(self,r_pnt):
        return self.parent_elem.get_dNr_mtx(r_pnt)
    
    def get_dNr_geo_mtx(self, r_pnt):
        return self.parent_elem.get_dNr_geo_mtx(r_pnt)
    
    def get_dNr_mtx(self,r_pnt):
        return self.parent_elem.get_dNr_mtx(r_pnt)
    
    def get_N_geo_mtx(self, r_pnt):
        return self.parent_elem.get_N_geo_mtx(r_pnt)