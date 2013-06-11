'''
Created on Jun 29, 2010

@author: alexander
'''
#import itertools
from itertools import product, chain
#itertools


gamma_list = [ [1, 2], [2, 3], [3, 4] ]


gamma_tuple = [1, 2], [2, 3], [3, 4]
gamma_test = chain( gamma_tuple )

print 'gamma_list', gamma_list
print 'gamma_tuple', gamma_tuple


print 'XXX1', list( product( ( [1, 2], [2, 3], [3, 4] ) ) )

print 'XXX2', list( product( gamma_list ) )

print 'XXX3', list( product( gamma_tuple ) )

print 'XXX4', list( product( gamma_test ) )

#
#def product( *args, **kwds ):
#
#    # product('ABCD', 'xy') --> Ax Ay Bx By Cx Cy Dx Dy
#    # product(range(2), repeat=3) --> 000 001 010 011 100 101 110 111
#
#    pools = map( tuple, args ) * kwds.get( 'repeat', 1 )
#    print 'args', args
##
##    pools = map( tuple, gamma_list )
#
#    print 'pools', pools
#
#    result = [[]]
#    for pool in pools:
#        print 'pool', pool
#        result = [x + [y] for x in result for y in pool]
#        print 'result', result
#
#    for prod in result:
#        yield tuple( prod )
#
#
#list( product( gamma_list ) )
#
#
#
#
print list( chain( gamma_list ) )
