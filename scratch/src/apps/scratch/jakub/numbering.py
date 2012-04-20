from numpy import arange

a = arange(0,100,18)
b = arange(16)
c = a[:,None] + b
print c.flatten()
