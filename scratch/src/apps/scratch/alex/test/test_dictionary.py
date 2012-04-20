from numpy import *



etype_list = ['G', 'E']
vcomp_list = ['N', 'T']

dictionary = { 'G_mac_N' : 1,  'E_mac_N' : 2, 'G_mac_T' : 3, 'E_mac_T' : 4 }


result_list_mac = {}
    
for etype in etype_list:
    for vcomp in vcomp_list:
    
        rtrace_key = etype + '_mac_' + vcomp
        print 'rtrace_key', rtrace_key
        
        mp_se_ee = dictionary[ rtrace_key ]
        xdata = copy( mp_se_ee )
        result_list_mac[ rtrace_key ] = mp_se_ee

print 'result_list_mac', result_list_mac

print result_list_mac.keys()
print result_list_mac.values()


a = "lkj"
print a
print type(a)

b = array(result_list_mac.keys())
print type(b)
print b

#print 'E_t_mic = %.10f' %(E_t )
