'''
Created on Jul 1, 2010

@author: alexander
'''
d = {'A':1, 'B':2}
c = {'C':22}
for k, v in d.iteritems():
    print k
    print v
    c[k] = v

print c
