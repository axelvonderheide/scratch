
from scipy.linalg import *
from numpy import *
from time import *



n_time = 1
i = 0
diff_t2_t1 = 0
diff_t4_t3 = 0

while i < n_time:
     
    psi = array([[1,2],[3, 4]])
    delta = array([[1,2],[3, 4]])
    M = zeros((2,2,2,2), dtype = float)
    
    
    t_loop = 0.
    for i in range(n_time):
    
        t1 = time()
        for i in range(0,2):
            for j in range(0,2):
                for k in range(0,2):
                    for l in range(0,2):
                        M[i,j,k,l] = 0.25 * ( dot(psi[i,k],delta[j,l]) + dot(psi[i,l],delta[j,k]) + dot(psi[j,k],delta[i,l]) + dot(psi[j,l],delta[i,k]))
        
        t2 = time()
        diff_t2_t1 = diff_t2_t1 + t2 - t1
        
        
#        print "time1 = ", t1
#        print "time2 = ", t2
#        print "t2 - t1 = ", t2 - t1
#        print "                   "
        
        
        t3 = time()
        
        M_ijkl = outer(psi, delta).reshape(2,2,2,2)
        M_ikjl = M_ijkl.swapaxes(1,2)

        M_iljk = M_ikjl.swapaxes(2,3)
        
        M_jkil = M_iljk.swapaxes(1,2)
        
        M_jlik = M_jkil.swapaxes(2,3)
        
        
        M_ijkl = 0.25 * ( M_ikjl + M_iljk + M_jkil + M_jlik )
        M_ijkl_arr = M_ijkl.reshape(2,2,2,2)
        
        print 'outer::::::::::::::::::::::: \n'
        print 'M_ikjl \n'
        print  M_ikjl 
        print 'M_iljk \n'
        print  M_iljk 
        print 'M_jkil \n'
        print  M_jkil
        print 'M_jlik \n'
        print  M_jlik
        
        
        t4 = time()

        diff_t4_t3 = diff_t4_t3 + t4 - t3

#        print "time3 = ", t3
#        print "time4 = ", t4
#        print "t4 - t3 = ", t4 - t3

        i = i + 1

        print "diff_t2_t1",diff_t2_t1/i
        print "diff_t4_t3",diff_t4_t3/i

        a = outer(psi, delta).reshape(2,2,2,2).swapaxes(1,2) 
        b = ((outer(psi,delta).reshape(2,2,2,2)).swapaxes(1,2)).swapaxes(2,3)  
        c = (outer(psi, delta).reshape(2,2,2,2).swapaxes(0,2)).swapaxes(1,2) 
        d = (((outer(psi, delta).reshape(2,2,2,2)).swapaxes(0,2)).swapaxes(1,2)).swapaxes(2,3) 

        print 'a \n', a-M_ikjl
        print 'b \n', b-M_iljk
        print 'c \n', c-M_jkil
        print 'd \n', d-M_jlik


        #y = 0.25 * ( outer(psi, delta) + outer(psi.transpose(),delta.transpose()) + outer(psi.transpose(), delta) + outer(psi, delta.transpose()) )
        #zz = y.reshape(2,2,2,2)


#        y = 0.25 * (( outer(psi, delta).reshape(2,2,2,2)).swapaxes(1,2) + 
#                    ((outer(psi,delta).reshape(2,2,2,2)).swapaxes(1,2)).swapaxes(2,3) + 
#                    ((outer(psi, delta).reshape(2,2,2,2)).swapaxes(0,2)).swapaxes(1,2) + 
#                    (((outer(psi, delta).reshape(2,2,2,2)).swapaxes(0,2)).swapaxes(1,2)).swapaxes(2,3) )
#        z = y.reshape(2,2,2,2)


#print "diff_t2_t1",diff_t2_t1/n_time
#print "diff_t4_t3",diff_t4_t3/n_time


