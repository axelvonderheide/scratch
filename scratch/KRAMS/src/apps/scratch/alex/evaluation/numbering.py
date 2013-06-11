from numpy import arange

#a = arange(1,22,17)
n = 3
a = arange(1,(n-1)*17+2,17)
n = 3
a = arange(1,n*17,17)
print "\na ",a

b = arange(16)
print "\nb ",b

c = a[:,None] + b
print "\nc ",c


c.flatten()


#aa = arange(1,10+17,17)
#bb = arange(16)
#cc = aa[:,None] + bb
#dd = cc.flatten()
#row_no = dd[:]
#print 'NUMBER2', row_no.shape[0]
