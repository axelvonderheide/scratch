
from numpy import allclose, arange, eye, linalg, ones, ix_, array, zeros
from scipy import linsolve, sparse,linalg

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
        
if __name__ == '__main__':
    mtx = SysDenseMtx( 10 )
    #mtx[2,3] = 10
    #mtx[:,8] = -10
    #print mtx[2,3]
    #print mtx
    #submtx = ones( (2,2), dtype = 'float_' )
    submtx = array([[1.,-1.],[-1.,1.]])
    print submtx
    for i in range(9):
        mtx[ ix_((i,i+1),(i,i+1) ) ] += submtx
    
    #print mtx[:,8]
    rhs = zeros(10)
    rhs[9]= 1.
    mtx[0,:] = 0.
    mtx[:,0] = 0.
    mtx[0,0] = 1.
    print mtx
    print mtx.solve(rhs)
    
    
    
    
    