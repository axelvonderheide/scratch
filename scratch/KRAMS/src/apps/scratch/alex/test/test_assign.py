from numpy import  *

e_list = [1,2,3]
print 'e_list', e_list

print 'len(e_list)', len(e_list)

a = zeros(len(e_list))
print 'a', a

a[:]=e_list
print 'a', a

b = array([7,7,7])
print 'b', b

c = max(a,b)
print 'c', c
