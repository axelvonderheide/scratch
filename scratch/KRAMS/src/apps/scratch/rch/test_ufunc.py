

from numpy import array,frompyfunc

arr = array([[[1,2,3],[2,3,4]],
             [[2,3,4],[3,4,2]]])

def getit( g ):
    return g*2, g*3

print map( getit, arr )