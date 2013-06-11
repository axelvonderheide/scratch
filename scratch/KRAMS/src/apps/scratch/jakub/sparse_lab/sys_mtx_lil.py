
from numpy import allclose, arange, eye, linalg, ones, ix_, array, zeros
from scipy import linsolve, sparse,linalg

from numpy import zeros, isscalar, real, imag, asarray, asmatrix, matrix, \
                  ndarray, amax, amin, rank, conj, searchsorted, ndarray,   \
                  less, where, greater, array, transpose, empty, ones, \
                  arange, shape, intc, clip, prod, unravel_index

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
        return linalg.solve( self.mtx, rhs )
    
class SysSparseMtx( object ):
    
    def __init__(self, ndofs ):
        self.mtx = sparse.lil_matrix((ndofs,ndofs))
        
    def __setitem__(self, idx, value ):
        self.mtx[idx] = value

    def __getitem__(self, idx):
        return self.mtx[idx]
    
    def __str__(self):
        return str( self.mtx )
    
    def solve(self, rhs):
        return linsolve.spsolve( self.mtx, rhs )   

        
if __name__ == '__main__':
    ndofs = 10
    submtx = array([[1.,-1.],[-1.,1.]])
    rhs = ones(ndofs)
    #rhs = zeros(ndofs)
    #rhs[-1] = 1.
    ix_list = [] #create index array
    for i in range(ndofs-1):
        ix_list.append([i, i+1])

#    tf_begin = time()
#    mtx_full = SysDenseMtx(ndofs)
#    for ix in ix_list:
#        mtx_full[ix_(ix,ix)] += submtx
# 
#    mtx_full[0,:] = 0.
#    mtx_full[:,0] = 0.
#    mtx_full[0,0] = 1.
#    #print "U_vct_full ",mtx_full.solve(rhs)
#    mtx_full.solve(rhs)
#    tf_end = time()
#    diff_f = tf_end -tf_begin
#    print "Full Matrix: %8.2f sec" % diff_f 
 
    ts_begin = time()
    mtx_sps = SysSparseMtx(ndofs)
    for ix in ix_list:
        mtx_temp = sparse.lil_matrix((ndofs,ndofs))
        for a in range(submtx.shape[0]):
            mtx_temp.rows[ix[a]]= ix
            mtx_temp.data[ix[a]]= submtx[a]
        mtx_sps.mtx = mtx_sps.mtx +  mtx_temp
 
    mtx_sps[0,:] = 0.
    mtx_sps[:,0] = 0.
    mtx_sps[0,0] = 1.
    #print "U_vct_sps ",mtx_sps.solve(rhs)
    mtx_sps.solve(rhs)
    ts_end = time()
    diff_s = ts_end -ts_begin
    print "Sparse Matrix: %8.2f sec" % diff_s 
    
    
    
    