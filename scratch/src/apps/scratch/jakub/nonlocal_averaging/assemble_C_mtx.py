'''
Created on Apr 7, 2010

@author: jakub
'''

from enthought.traits.api import \
     Array, Bool, Callable, Enum, Float, HasTraits, Interface, implements, \
     Int, Trait, Str, Enum, Callable, List, TraitDict, Any, \
     on_trait_change, Tuple, WeakRef, Delegate, Property, cached_property,\
     Instance, Dict
     
class AveragingFunction(HasTraits):
    radius = Float(.2)
    correction = Bool(True) 
    
    def set_center_old(self,X):
        self.center = X
    
    def get_value(self,dist):
        raise NotImplementedError
        
class QuarticAF(AveragingFunction):
    
    def get_value(self, dist):
        if dist > self.radius:
            return 0.
        else:
            return ((1-(dist/self.radius)**2)**2)**2 
    
    def get_value_old(self,X):
        value = (1-(norm(self.center-X)/self.radius)**2)**2 
        if value > 0:
            return value**2
        else:
            return 0
        
class LinearAF(AveragingFunction):
    
    def get_value(self, dist):
        if dist > self.radius:
            return 0.
        else:
            return (1./self.radius) * (1-dist/self.radius)
        
    def get_value_old(self,X):
        value = (1./self.radius) * (1-norm(self.center-X)/self.radius)
        if value > 0.:
            return value
        else:
            return 0
     
################################
# copy of the C_mtx property
################################
        
    C_mtx = Property
    @cached_property
    def _get_C_mtx(self):
        #--------------
        #Averaging
        #------------------
        t1 = time()
        n_dofs = self.sdomain.n_dofs
        a_fn = frompyfunc(self.averaging.get_value, 1, 1)
        dim = self.fets_eval.n_nodal_dofs #TODO:find better way
        #TODO:hack, works just with one fe_grid, 
        #make it work for random number of fe_grids
        point_grid = self.sdomain.fe_subgrids[0].dof_grid.cell_grid.point_grid
        #TODO:check if the array is not directly aviable in cell_grid
        X_pnt = point_grid.reshape((dim,-1)).T   
        ip_coords = self.fets_eval.ip_coords
        n_ip = len(ip_coords)
        ip_weights = self.fets_eval.ip_weights
        data, row, col = [],[],[]  
        #if the same integration scheme is used the N matrices can be evaluated 
        #just once - shape(n_ip,dims, dofs)
        e_ip_N_mtx = array([self.fets_eval.get_N_mtx(r_pnt) for r_pnt in ip_coords])
        n_elem_dofs = e_ip_N_mtx.shape[2]
        #Loop over all nodes 
        for i,x_pnt in enumerate(X_pnt):
            #print 'x_pnt ',x_pnt
            #Use LS to find the elements inside the radius
            #TODO: make it for all subgrids
            #TODO:check leaking
            level_set = self.sdomain.fe_subgrids[0]\
                        ['(X-%(x)f)**2 + (Y-%(y)f)**2 - %(r)f**2'\
                        % {'x': x_pnt[0], 'y': x_pnt[1],'r':self.averaging.radius}]
            active_elems = hstack((level_set.elems,level_set.neg_elems))
            n_a_elems = active_elems.size
            if n_a_elems == 0:
                'WARNING - Radius too small!'
            #Generate the IP coords for all active elems            
            #element transformation have to be used due distorted meshes
            elems_ip_coords = zeros((n_a_elems,n_ip, dim), dtype = float)
            elems_dof_map = zeros((n_a_elems,n_elem_dofs), dtype = int)
            for j,e_id in enumerate(active_elems):
                X_mtx = self.sdomain.fe_subgrids[0].elements[e_id].points
                elems_ip_coords[j,:] = \
                            self.fets_eval.get_vtk_r_glb_arr(X_mtx, ip_coords)
                elems_dof_map[j,:] = \
                        self.sdomain.fe_subgrids[0].elements[e_id].dofs
            arm = elems_ip_coords - x_pnt[None,None,:]
            dist = cdist(x_pnt[None,:], elems_ip_coords.\
                         reshape(n_a_elems*n_ip,dim)).\
                         reshape(n_a_elems,n_ip)
            alpha = a_fn(dist)#(elems,n_ip)
            J_det = self.J_det_grid[active_elems]#(elems, n_ip)
            values = (ip_weights.T*J_det) * alpha#(elems, n_ip)
            r_00 = values.sum()
            if self.averaging.correction:
                r_01 = np_sum(np_sum(values[...,None]*arm, axis = 1),axis = 0)
                outer_p = arm[...,None,:]*arm[...,None]#outer product (elems,n_ip, dims, dims)
                #TODO:is summing over more dims ad once posible?
                R_11 = np_sum(np_sum(values[...,None,None]* outer_p,axis = 1),axis = 0)
                #evaluate the correcting factors
                A = zeros((dim+1,dim+1))
                A[0,0] = r_00
                A[0,1:] = r_01
                A[1:,0] = r_01
                A[1:,1:] = R_11
                b = zeros((dim+1), dtype = float )
                b[0] = 1.
                params = linalg.solve(A,b)
                #print 'check ',dot(A,params)
                p0 = params[0]
                p1 = params[1:]
                m_correction = np_sum(p1[None,None,:]*arm, axis = 2)#(elem,n_ip)
                c_elem = np_sum((values*(p0 + m_correction))[...,None,None]*\
                         e_ip_N_mtx[None,...],axis =1)
            else:
                p0 = 1./r_00
                c_elem = np_sum((values * p0)[...,None,None]*e_ip_N_mtx[None,...],axis = 1)

            data.append(c_elem.flatten())
            col.append(tile(elems_dof_map,dim).flatten())
            row.append(((i*dim)+zeros_like(c_elem)+\
                        arange(dim, dtype = int)[None,:,None]).flatten())
        data = hstack(data)
        row =  hstack(row)
        col = hstack(col)

        t2 = time()
        diff = t2 - t1
        print "Averaging Matrix: %8.2f sec" % diff 
        return coo_matrix((data,(row,col)),shape = (n_dofs, n_dofs), dtype = float_ ).tocsr()

################################
# multiplication in corr_pred
################################

        u_avg = C_mtx * u.T
