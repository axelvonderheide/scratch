from numpy import *

aa = arange(3)

di = { aa:10, "b":20, "c":30}

#di = { "a":10, "b":20, "c":30}
print 'di' , di

aa = arange(3)
print 'aa' ,aa

print 'aa' , list(aa)


print 'di.keys()' , di.keys()
print 'di.values()' , di.values()

#di.values = list(aa)