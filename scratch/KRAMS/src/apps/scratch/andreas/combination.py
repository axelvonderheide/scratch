'''
Created on Sep 15, 2009

@author: Andreas
'''
from enthought.traits.api import \
    HasTraits, Directory, List, Int, Float, Any, Enum, \
    on_trait_change, File, Constant, Instance, Trait, \
    Array, Str, Property, cached_property, WeakRef, \
    Dict, Button, Color

from numpy import \
    copy, ones, arange, vstack,hstack, shape, array, append, transpose, \
    arange, reshape, c_, newaxis, sum, insert, min, max, argmin, argmax, sort

class combination( HasTraits ):
    
    quantity_g = 1
    quantity_q = 4
    
    # psi values, 0, 1, 2
        
    psi_q_1 = array([1.0,0.4,0.0])
    psi_q_2 = array([0.8,0.4,0.0])
    psi_q_3 = array([0.7,0.4,0.0])
    psi_q_4 = array([0.6,0.5,0.0])

    psi_array= vstack((psi_q_1,psi_q_2,psi_q_3,psi_q_4)) 

    def get_comb_uls(self):
        """possible combination for uls"""
        
        # gamma values
        
        gamma_list = [[1.35,1.0],[1.5,0.0],[1.5,0.0],[1.5,0.0],[1.5,0.0]]

        
        # combination of all possibilities
        def product(*args, **kwds):
            # product('ABCD', 'xy') --> Ax Ay Bx By Cx Cy Dx Dy
            # product(range(2), repeat=3) --> 000 001 010 011 100 101 110 111
            pools = map(tuple, args) * kwds.get('repeat', 1)
            result = [[]]
            for pool in pools:
                result = [x+[y] for x in result for y in pool]
            for prod in result:
                yield tuple(prod)


        r = perm(gamma_list)

        l=[]
        for line in r:
            list_str = line.split(";")[:-1]
            l.append(list_str)

        comb = array(l, dtype = 'float_')
        
        # load case combinations 
        # for different leading load cases including psi, values
        
        lead_1 = copy(self.psi_array[:,0])
        lead_1[0] = 1.0
    
        lead_2 = copy(self.psi_array[:,0])
        lead_2[1] = 1.0
        
        lead_3 = copy(self.psi_array[:,0])
        lead_3[2] = 1.0
        
        lead_4 = copy(self.psi_array[:,0])
        lead_4[3] = 1.0
        
        lead_1 = hstack((ones(self.quantity_g),lead_1))
        lead_2 = hstack((ones(self.quantity_g),lead_2))
        lead_3 = hstack((ones(self.quantity_g),lead_3))
        lead_4 = hstack((ones(self.quantity_g),lead_4))    
        
        comb_1 = comb * lead_1
        comb_2 = comb * lead_2
        comb_3 = comb * lead_3
        comb_4 = comb * lead_4
        
        comb = vstack((comb_1,comb_2,comb_3,comb_4))
        
        multiple_i =array([])
        multiple_j =array([])
        
        for i in range(0,len(comb)):
            for j in range(0,len(comb)):
                if i!=j and (comb[i] == comb[j]).all():
                    multiple_i = append(multiple_i,i)
                    multiple_j = append(multiple_j,j)
        print sort(multiple_j)
        print sort(multiple_i)
        return comb
    
    def get_comb_sls(self):
        """possible combination for sls"""
        
        # gamma values
        
        gamma_g_1 = array([1.0,1.0])
        gamma_q_1 = array ([1.0,0.0])
        gamma_q_2 = array ([1.0,0.0])
        gamma_q_3 = array ([1.0,0.0])
        gamma_q_4 = array ([1.0,0.0])
        
        
        # combination of all possibilities
        
        comb = array    
        n=0.0                
        for y_g in gamma_g_1:
            for y_q_1 in gamma_q_1:
                for y_q_2 in gamma_q_2:
                    for y_q_3 in gamma_q_3:
                        for y_q_4 in gamma_q_4:
                            a=array([[y_g,y_q_1,y_q_2,y_q_3,y_q_4]])
                            if n == 0:
                                comb = a
                            else:
                                comb = append(comb,a, axis =0)
                            n = n+1.0

        # load case combinations 
        # for different leading load cases including psi, values

        lead_1 = copy(self.psi_array[:,2])
        lead_1[0] = 1.0
        
        lead_2 = copy(self.psi_array[:,2])
        lead_2[1] = 1.0
        
        lead_3 = copy(self.psi_array[:,2])
        lead_3[2] = 1.0
        
        lead_4 = copy(self.psi_array[:,2])
        lead_4[3] = 1.0
        
        lead_1 = hstack((ones(self.quantity_g),lead_1))
        lead_2 = hstack((ones(self.quantity_g),lead_2))
        lead_3 = hstack((ones(self.quantity_g),lead_3))
        lead_4 = hstack((ones(self.quantity_g),lead_4))    
        
        comb_1 = comb * lead_1
        comb_2 = comb * lead_2
        comb_3 = comb * lead_3
        comb_4 = comb * lead_4
        
        
        print n
        comb = vstack((comb_1,comb_2,comb_3,comb_4))
        
        return comb
    
    def get_lcc(self,lc, combin):
        """
        estimation of load case combination
        given load cases need to have the following shape:
        (# load cases, # elements, # stress resultants)
        f.e. lc =array([load1, load2, load3, load4, load5])
        """
        #combine all load cases
        #
        
        lcc_arr = sum( lc  * (combin[:,:,newaxis,newaxis]),axis=1.0)
        
        # create an array for identification of the load case combination
        #
        
        n_lcc, n_elem = shape(lcc_arr)[0:2]
        
        idx = transpose(ones((n_elem,n_lcc))* (arange(n_lcc)))[:,:,newaxis]
        
        # stack index into lcc as last value 
        # after internal forces (mx, my, mxy,...., idx)
        #
        
        lcc_i  = transpose(vstack((transpose(lcc_arr), transpose(idx))))
        
        print "lcc",shape(lcc_arr)
        print "idx",shape(idx)
        print "lcc_i",shape(lcc_i)
        
        # final shape of lcc_i 
        # ( n_lc , n_elem, n_internal forces+idx(len=1) )
        #
        
        return lcc_i
    
    def get_min_max(self,lcc_arr_i):
        """
        for traceability of combination
        """
        
        arg_min     =   argmin(lcc_arr_i[:,:,:-1],axis=0)
        print " arg_min", shape(arg_min), arg_min
        min_val     =   min(lcc_arr_i[:,:,:-1],axis=0)
        arg_max     =   argmax(lcc_arr_i[:,:,:-1],axis=0)   
        max_val     =   max(lcc_arr_i[:,:,:-1],axis=0)
        
        n_elemens = shape(lcc_arr_i)[1]
        
        min_max_arr = vstack((min_val, arg_min, max_val, arg_max))
        min_max_arr = min_max_arr.reshape(shape(min_max_arr)[0]
                                          /n_elemens,n_elemens,6) 
        
        return min_max_arr
        
if __name__ == '__main__':
    com = combination()
    
    com_uls = com.get_comb_uls()

    com_sls = com.get_comb_sls()
    
    # example data
    #
    
    data1 = arange(12).reshape(2,6)+1.0
    data2 = arange(12).reshape(2,6)*0
    data3 = arange(12).reshape(2,6)*0
    data4 = arange(12).reshape(2,6)*0
    data5 = arange(12).reshape(2,6)*0
    
    # array of all load cases for method get_lcc
    #
    
    #lc_arr =array([data1, data2, data3, data4, data5])
    
    
    #lcc_arr_i = com.get_lcc(lc_arr, com_uls)
    
    #min_max = com.get_min_max(lcc_arr_i)
    #print shape(min_max)
    #min_mx = argmin(lcc_i[:,:,:],axis=0)
    #print min_mx
    

    