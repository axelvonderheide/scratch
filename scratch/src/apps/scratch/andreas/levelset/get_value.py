'''
Created on Sep 4, 2009

@author: Andreas
'''
#distance estimation using euclidean definition for mesh points using python2.6

from scipy.spatial.distance import cdist
from numpy import min       #A:Isnt this part of the main program anyway??
                            #J.this overloads the default min, this one is 
                            #faster and more powerful on arrays
#from scipy import transpose#J:this one is included in numpy arrays, you can just go A.T
#import numarray             #A: didnt manage to import only the function numarray.minimum.reduce


def get_value(imput_pts,grid):
    #resultall=result=[]#J: you don't need to initiate in python, unless you just wanna append later on
                        #intiating two params like that is dangerous, because they stay linked, 
                        #so if you change one of them, the other changes too
    resultall=cdist(imput_pts,grid)
    #result=numarray.minimum.reduce(resultall.T)    #A:Transpose - Inefficient?
                                                    #j:transpose has no penalty, it just reads the array
                                                    #in different order
    result = min(resultall, axis=0)#J: this seems to be working
    return result

#test
if __name__ == '__main__':
    print 'distance ',get_value([[0.,0.],[0,2.]],[[0.,0.5],[0,1.],[0,1.5]])

#sugestions:
#min() does that all, the numpy functions obviously do have loops inside, but written in fortran
#http://mail.python.org/pipermail/python-list/2004-March/251218.html

# answer:
# because cdist gives me a matrix or vector i need to solve for the maximum of this vector or matrix
# to avoid loops i inplemented numarray.minimum.reduce, what if this function is using a loop in itself???


#future steps:
#1. create object SifnDistFn with param input_pts and your method
#2. use numpy random module and mgrid functions, to generate 
#iregular points / regular grid to have some bigger samples to test time on

# Nice weekend to you too
# Greetings Jakub