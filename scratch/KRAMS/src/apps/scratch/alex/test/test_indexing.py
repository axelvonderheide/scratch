from numpy import *
#
#a = array([10, 2, -4, 5])
#print 'a', a
#
#a_ix = argsort(a)
#print 'a_ix', a_ix
#
#a_sorted = a[a_ix]
#print 'a_sorted', a_sorted
#
#print 'XXX', a[::-1]
#
#a_list = list(a_sorted)
#a_list.reverse()
#print 'a_', a_list

a = arange(100)
for i in range(100):
    print 1 + i * 16
#    a[i]= 1 + i*16
#    a[i]= 2 + i*16
#    a[i]= 3 + i*16
#    a[i]= 4 + i*16
#    a[i]= 5 + i*16
#    a[i]= 6 + i*16
#    a[i]= 7 + i*16
#    a[i]= 8 + i*16
#    a[i]= 9 + i*16
#    a[i]= 10 + i*16
#    a[i]= 11 + i*16
#    a[i]= 12 + i*16
#    a[i]= 13 + i*16
#    a[i]= 14 + i*16
#    a[i]= 15 + i*16
#    a[i]= 16 + i*16

mx = 0.56
m2 = -0.01
mxy = -0.02

print arctan(mxy/(m2-mx))
