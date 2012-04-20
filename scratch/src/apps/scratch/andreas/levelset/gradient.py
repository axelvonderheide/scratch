'''
Created on Sep 16, 2009

@author: jakub
'''
from sympy import diff,symbols

x,y = symbols('xy')


f = x**2+y**2 - 25.

df = diff(f,x,1)

f1 = f/df

print 'df ',df
print 'f1 ',f1