
from numpy import allclose, arange, eye, linalg, ones, ix_, array, zeros, \
                hstack,meshgrid, vstack
from scipy import linsolve, sparse,linalg
from time import time

class SysDenseMtx( object ):
    
    def __init__(self, ndofs ):
        self.mtx = zeros( (ndofs,ndofs), dtype = 'float_' )
        
    def __setitem__(self, idx, value ):
        self.mtx[idx] = value

    def __getitem__(self, idx):
        return self.mtx[idx]
    
    def __str__(self):
        return str( self.mtx )

    def solve(self, rhs):
        tf_solve_s = time()
        u_vct = linalg.solve( self.mtx, rhs )
        tf_solve_e = time()
        dif_solve = tf_solve_e - tf_solve_s
        print "Full Solve: %8.2f sec" %dif_solve
        return u_vct
    
class SysSparseMtx( object ):
    
    def __init__(self, n_dofs, data_l, x_l, y_l ):
        #self.mtx = sparse.sparse.coo_matrix((data,ij))
        self.n_dofs = n_dofs
        self.data_l = data_l
        self.x_l = x_l
        self.y_l = y_l
        
    
#    def set_nond_row_zero(self, a):
#        for i in range( len(self.data_l)-1, -1, -1):# loop has to run backwards
#            if self.x_l[i] == a and self.y_l[i] != a:
#                self.x_l.pop(i)
#                self.y_l.pop(i)
#                self.data_l.pop(i)
#                
#    def set_nond_col_zero(self, a):
#        for i in range(len(self.data_l)-1, -1, -1):# loop has to run backwards
#            if self.y_l[i] == a and self.x_l[i] != a:
#                self.x_l.pop(i)
#                self.y_l.pop(i)
#                self.data_l.pop(i)
                
    def set_nond_zero(self, a):
        for i in range(len(self.data_l)-1, -1, -1):# loop has to run backwards
            if self.x_l[i] == a or self.y_l[i] == a:
                if self.y_l[i] == a and self.x_l[i] == a:
                    self.data_l[i] = - self.data_l[i]
                else:
                    self.x_l.pop(i)
                    self.y_l.pop(i)
                    self.data_l.pop(i)
    
    def get_column(self, a):                
        K_ac = zeros(self.n_dofs)
        for i in range(len(self.data_l)):
            if self.y_l[i] == a:
                K_ac[self.x_l[i]] += self.data_l[i]
        return K_ac
    
    def get_diag_elem(self, a):
        K_aa = 0.
        for i in range(len(self.data_l)):
            if self.x_l[i] == a and  self.y_l[i] == a:
                K_aa += self.data_l[i]
        return K_aa
        
    def solve(self, rhs):
        ij = vstack((self.x_l,self.y_l))
        #print "ij", ij
        #print "data ", self.data_l
        mtx = sparse.coo_matrix((self.data_l,ij))
        mtx.tocsr()
        ts_solve_s = time()
        u_vct = linsolve.spsolve( mtx, rhs )
        ts_solve_e = time()
        dif_solve = ts_solve_e - ts_solve_s
        print "Sparse Solve: %8.2f sec" %dif_solve
        return u_vct  

        
if __name__ == '__main__':
    n_dofs = 5000
    submtx = array([[1.,-1.],[-1.,1.]])
    #print "submtx",submtx
    ix_list = []
    for i in range(n_dofs-1):#generate indes matrix
        ix_list.append(array([i, i+1]))
    rhs = zeros(n_dofs)
    rhs[-1]= 1.
    
    tf_start = time()
    mtx = SysDenseMtx( n_dofs )
    for ix in ix_list:    
        mtx[ ix_(ix,ix)  ] += submtx
    
   

    mtx[0,:] = 0.
    mtx[:,0] = 0.
    mtx[0,0] = 1.
    #print "mtx ",mtx
    U_f=mtx.solve(rhs)
    print "U_f ",U_f
    tf_end = time()
    diff = tf_end - tf_start
    print "Full Matrix: %8.2f sec" %diff
    
    ts_start = time()
    x_l = []
    y_l = []
    data_l = []
    for ix in ix_list:
        X,Y = meshgrid(ix,ix)
        x_l.extend(X.flatten())
        y_l.extend(Y.flatten())
        data_l.extend(submtx.flatten())
    
    #print "X,Y ", x_l,y_l
    #print "data ", data_l
    #print "X,Y ",x_a,y_a
    #x_a = array([x_l]).flatten()
    #y_a = array([y_l]).flatten()
    #ix_arr = vstack((x_a,y_a))
    #print "ix_arr ", ix_arr
    #data = array(data_l).flatten()
    #print data.shape[0]
    sp_mtx = SysSparseMtx(n_dofs, data_l, x_l, y_l)
    #sp_mtx.set_nond_row_zero(0)
    #sp_mtx.set_nond_col_zero(0)
    sp_mtx.set_nond_zero(0)
    #K_aa =  sp_mtx[0,:]       
    #sp_mtx[0] = 0.
    #sp_mtx[:,0] = 0.
    #sp_mtx[0,0] = 1.
    #print "col ",sp_mtx.get_column(1)
    #print "K_aa ",sp_mtx.get_diag_elem(2)
    U_s = sp_mtx.solve(rhs)
    print "U_s ", U_s
    ts_end = time()
    difs = ts_end - ts_start
    print "Sparse Matrix: %8.2f sec" %difs

#    for i in range(5,0,-1):
#        print "i ",i