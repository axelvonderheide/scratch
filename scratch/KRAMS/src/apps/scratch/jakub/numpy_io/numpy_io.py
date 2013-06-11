from numpy import array, save, load
save('test', array([1,2,3]))
print load('test.npy')
print "done"