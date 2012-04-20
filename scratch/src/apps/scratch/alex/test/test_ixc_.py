from numpy import  ix_, setdiff1d, arange

'''
definition of a complementary index function extracting the complementary indices
and extracting the corresponding rows and colomns from the array a
'''

a = arange(24).reshape(4,6)

print 'a = \n', a
print '\n'

print 'setdiff1d(range(a.shape[0]),[2,5]) = \n', setdiff1d(range(a.shape[0]),[2,5])
print '\n'

print 'ix_([1,3],[2,5]) = \n', ix_([1,3],[2,5])
print '\n'
    
b = a[ix_([1,3],[2,5])]

bc = a[ix_(setdiff1d(range(a.shape[0]),[1,3]), \
           setdiff1d(range(a.shape[1]),[2,5]))]

bcc = a[ix_([0, 2],[0,1,3,4])]


print '\n'
print 'b', b

print '\n'
print 'bc', bc
print '\n'

print '\n'
print 'bcc', bcc
print '\n'

     
def ixc_(a_array, ix_list):
    ixc1 = setdiff1d( range(a_array.shape[0]), ix_list[0] )
    ixc2 = setdiff1d (range(a_array.shape[1]), ix_list[1] )
    ac = a[ix_(ixc1, ixc2)]
    return ac

print 'ixc_(a,[[1,3],[2,5]]) \n', ixc_(a, [[1,3],[2,5]])

a[0:self.shape[0]]
