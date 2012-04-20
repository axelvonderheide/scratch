from numpy import \
    shape, vstack, array, append, dtype,hstack

def product(*args):
    # product('ABCD', 'xy') --> Ax Ay Bx By Cx Cy Dx Dy
    # product(range(2), repeat=3) --> 000 001 010 011 100 101 110 111
    pools = map(tuple, args)
    result = [[]]
    for pool in pools:
        #print "pool", pool
        result = [x+[y] for x in result for y in pool]
        #print "result", result
    #print result

def product2(*args):
    # product('ABCD', 'xy') --> Ax Ay Bx By Cx Cy Dx Dy
    # product(range(2), repeat=3) --> 000 001 010 011 100 101 110 111
    pools = args
    print pools
    result = [[]]
    for pool in pools:
        print "pool", pool
        result = [x+[y] for x in result for y in pool]
        print "result", result
    print result
    
#r = product([1.35,1.0],[1.5,0],[1.5,0])
a= (1.35,1.0),(1.5,0),(1.5,0)

print a[:]
r2 = product2(([1.35,1.0],[1.5,0],[1.5,0]))



#print r
#print r2
#r = vstack(r)

#p#rint shape(r)