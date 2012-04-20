from numpy import  *

# numpy - broadcasting

#broadcasting b.shape to (2,2) (copy scalar to all)
a = array( [[1, 2], [3, 4]] )
print 'a \n', a
print 'a.shape ', a.shape
b = array( [1] )
print 'b \n', b
print 'b.shape', b.shape
a + b
print 'a+b \n', a + b
print ' \n'

#broadcasting b.shape to (2,2) (copy row to new rows)
b = array( [1, 2] )
print 'b \n', b
a + b
print 'a+b \n', a + b

#broadcasting a and b to shape (3,2) (works only with axis equal to 1!)
a = array( [[1, 2]] )
print 'a ', a
print 'a.shape ', a.shape
b = array( [[3], [4], [5]] )
print 'b.shape ', b.shape
print 'b ', b
a + b
print 'a+b \n', a + b
print ' \n'

# or use newaxis: this corresponds to adding "1," in shape
print 'newaxis:'
a = array( [1, 2] )
print 'a ', a
print 'a.shape ', a.shape
b = array( [3, 4, 5] )
print 'b.shape ', b.shape
print 'b ', b
c = a[newaxis, :] + b[:, newaxis]
print 'a+b ', c
print ' \n'


#
a = zeros( ( 3, 4, 5, 6, 7 ) )
print 'a.shape \n', a.shape
b = zeros( ( 7, ) )
print 'b.shape \n', b.shape
c = a + b
print 'c.shape', c.shape
print ' \n'

#
a = zeros( ( 3, 4, 5, 6 ) )
print 'a.shape \n', a.shape
b = zeros( ( 4, 6 ) )
print 'b.shape \n', b.shape
c = a + b[:, newaxis, :]
print 'c.shape', c.shape
print ' \n'

#
a = ones( ( 3, 4, 5, 6 ) )
print 'a.shape \n', a.shape
b = ones( ( 4, 6 ) )
print 'b.shape \n', b.shape
c = a + b[:, newaxis, :]
print 'c.shape', c.shape
print ' \n'
print 'c', c


a = arange( 6 ).reshape( 3, 2 )
print 'a.shape \n', a.shape
b = arange( 6 ).reshape( 3, 2 )
print 'b.shape \n', b.shape
c = a * b
print 'c.shape', c.shape
print ' \n'
print 'c', c


a = arange( 6 * 4 ).reshape( 3, 2, 4 )
print 'a.shape \n', a.shape
b = arange( 3 * 4 ).reshape( 3, 1, 4 )
print 'b.shape \n', b.shape
c = hstack( [ a, b ] )
print 'c.shape', c.shape
print ' \n'
print 'c', c

