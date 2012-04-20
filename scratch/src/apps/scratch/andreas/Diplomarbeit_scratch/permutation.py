from numpy import \
    shape, vstack, array, append, dtype,hstack

def product( args ):
    """
    Get all possible permutations of the security factors
    without changing the order of the loading cases.
    The method corresponds to the build-in function 'itertools.product'.
    Instead of returning a generator object a list of all 
    possible permutations is returned. As argument a list of list 
    needs to be defined. In the original version of 'itertools.product' 
    the function takes a tuple as argument ("*args"). 
    """
    pools = map( tuple, args ) #within original version args defined as *args
    result = [[]]
    for pool in pools:
        result = [x + [y] for x in result for y in pool]
    return result
    
    
#r = product([1.35,1.0],[1.5,0],[1.5,0])
a= (1.35,1.0),(1.5,0),(1.5,0),(1.5,1.0)
r2 = product(([1.35,1.0],[1.5,0],[1.5,0],[1.5,1.0]))

#print r
print len(r2)
#r = vstack(r)

#p#rint shape(r)