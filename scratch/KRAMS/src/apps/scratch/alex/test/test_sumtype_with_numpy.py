from scipy.linalg import *
from numpy import *
from time import *

'''check if with numpy functionality (i.e. using "outer()", "reshape", and "swapaxes()" )
   the same results are obtained then with the sum-type symmetrization formula stated in 
   Eq. [21] in "Comments on microplane theory" by Jirasek.
   The file shows that the same result is obtained and that a speed-up of factor 13 can be achieved!
''' 

### measure and compare the time needed for the two different approaches     
# time needed applying formula
t_numpy = 0  
# time needed applying numoy functionality
t_loop  = 0  


### arbitrary matrices for testing:
psi = array([[1,2],[3, 4]])
delta = array([[5,6],[7, 8]])

# 4th order tensors for components of the sum stated in the formula:
MM      = zeros((2,2,2,2), dtype = float)
MM_ikjl = zeros((2,2,2,2), dtype = float)
MM_iljk = zeros((2,2,2,2), dtype = float)
MM_jkil = zeros((2,2,2,2), dtype = float)
MM_jlik = zeros((2,2,2,2), dtype = float)
 
    
#for tt in range(10000):    
for tt in range(1):    
    
    # applying the formula Eq.[21]: 
    t1 = time()
    for i in range(0,2):
        for j in range(0,2):
            for k in range(0,2):
                for l in range(0,2):
                    MM[i,j,k,l] = 0.25 * ( dot(psi[i,k],delta[j,l]) + dot(psi[i,l],delta[j,k]) + \
                                           dot(psi[j,k],delta[i,l]) + dot(psi[j,l],delta[i,k]))
                    MM_ikjl[i,j,k,l] = psi[i,k]*delta[j,l]   
                    MM_iljk[i,j,k,l] = psi[i,l]*delta[j,k]
                    MM_jkil[i,j,k,l] = psi[j,k]*delta[i,l]
                    MM_jlik[i,j,k,l] = psi[j,l]*delta[i,k]

    t2 = time()
    t_loop += t2-t1

       
    # applying numpy functionality:
    t1 = time()
    M_ijkl = outer(psi, delta).reshape(2,2,2,2)
    Mc_ijkl = copy(M_ijkl)
        
    # it seams that the index notation always referse to the initial indexing (i=0,j=1,k=2,l=3)
    # Note: swapaxes returns a reference not a copy!!!
    M_ikjl = M_ijkl.swapaxes(1,2)
    M_iljk = M_ikjl.swapaxes(2,3)
    M_jlik = M_iljk.swapaxes(0,1)
    M_jkil = M_jlik.swapaxes(2,3)
#    # NOTE: this does NOT yield the same result 
#    M_ikjl = M_ijkl.swapaxes(1,2)
#    M_iljk = Mc_ikjl.swapaxes(1,3)
#    M_jlik = Mc_iljk.swapaxes(0,2)
#    M_jkil = Mc_jlik.swapaxes(1,3)

    # NOTE: this does yield the same result 
    Mc_ikjl_ = Mc_ijkl.swapaxes(1,2)
    Mc_ikjl = copy(Mc_ikjl_)
    Mc_iljk = Mc_ikjl.swapaxes(1,3)
    Mc_jlik = Mc_iljk.swapaxes(0,2)
    Mc_jkil = Mc_jlik.swapaxes(1,3)


    aaa = M_iljk == Mc_iljk
    print 'is this equal?', aaa.all()

#    M_sum = 0.25 * ( M_ikjl + M_iljk + M_jkil + M_jlik )
#
#    t2 = time()
#    t_numpy += t2-t1
#
#
## print and compare the results obtained with the different approaches.
## They both yield the same result!
#
#print 't_loop :', t_loop
#print 't_numpy:', t_numpy
## speed-up by using numpy-functionality instead of loop: about 13 times faster
#print 't_loop/t_numpy:', t_loop/t_numpy
#
#print 'formula::::::::::::::::::::::: \n'
#print 'MM_ikjl \n'
#print  MM_ikjl 
#print 'MM_iljk \n'
#print  MM_iljk 
#print 'MM_jkil \n'
#print  MM_jkil
#print 'MM_jlik \n'
#print  MM_jlik
#print 'MM \n'
#print  MM
#
#print 'using numpy::::::::::::::::::::::: \n'
#print 'M_ikjl \n'
#print  M_ikjl 
#print 'M_iljk \n'
#print  M_iljk 
#print 'M_jkil \n'
#print  M_jkil
#print 'M_jlik \n'
#print  M_jlik
#print 'M_sum \n'
#print  M_sum
#
#print 'difference::::::::::::::::::::::: \n'
#print 'Diff: M_ikjl \n'
#print  MM_ikjl - M_ikjl 
#print 'Diff: M_iljk \n'
#print  MM_iljk - M_iljk 
#print 'Diff: M_jkil \n'
#print  MM_jkil - M_jkil
#print 'Diff: M_jlik \n'
#print  MM_jlik - M_jlik
#print 'Diff: M_sum \n'
#print  M_sum - M_sum
#
#
#
#
#
#
