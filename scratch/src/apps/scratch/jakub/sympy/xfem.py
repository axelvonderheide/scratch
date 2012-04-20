'''
Created on Aug 10, 2009

@author: jakub
'''
from sympy import symbols, Matrix, Plot, diff, sign,sin,abs, integrate,pprint

#from sympy.functions import sign

from sympy.functions.special.delta_functions import Heaviside, DiracDelta


Xi_jump = 0.5
x,x_i,Xi = symbols('x,x_i,Xi') 

def SGN_S(x):
    return x/abs(x)

def SGN_H(x):
    return 2 * Heaviside(x) - 1.

def Jump_S(x,x_i,Xi) :
    return 1./2.*(SGN_S(x-Xi)-SGN_S(x_i-Xi))

def Jump_H(x,x_i,Xi):
    return 1./2.*(SGN_H(x-Xi)-SGN_H(x_i-Xi))

N_mtx_S = Matrix([1-x,x, (1-x)*Jump_S(x,0.,Xi_jump), x*Jump_S(x,1.,Xi_jump)])
N_mtx_H = Matrix([1-x,x, (1-x)*Jump_H(x,0.,Xi_jump), x*Jump_H(x,1.,Xi_jump)])

print "N_s\n",N_mtx_S
print "N_h\n",N_mtx_H

dx = lambda a: diff(a,x)
p = Plot()
for i,N in enumerate(N_mtx_S):
    p[i]=N,[x,0,1,100]
    
B_mtx_S = N_mtx_S.applyfunc(dx)
B_mtx_H = N_mtx_H.applyfunc(dx)

print 'B_s\n',B_mtx_S
print 'B_h\n',B_mtx_H

Ef, Em = symbols('Ef,Em') 
D_mtx = Matrix([Ef + Em *(1 - DiracDelta(Xi_jump))])

K_r_S = B_mtx_S * D_mtx*B_mtx_S.T
K_r_H = B_mtx_H * D_mtx*B_mtx_H.T

K_r_Sdata = K_r_S.subs(Ef,20).subs(Em,20)
K_r_Hdata = K_r_H.subs(Ef,20).subs(Em,20)

print 'K_r_s\n',K_r_S.shape
print 'K_r_h\n',K_r_H.shape



integr = lambda a:integrate(a,(x,0,1))
integr_l = lambda a:integrate(a,(x,0.,0.5))
integr_r = lambda a:integrate(a,(x,0.5,1.))

#print K_r_Sdata.applyfunc(integr_l)+K_r_Sdata.applyfunc(integr_r)
K_i_H = K_r_Hdata.applyfunc(integr)
print "KiH\n",K_i_H

K_i_Hc = K_r_Hdata.applyfunc(integr_l)+K_r_Hdata.applyfunc(integr_r)
print "KiHc\n",K_i_Hc

print "diff\n",K_i_H-K_i_Hc
#K_i_H.col_del(1)
#K_i_H.row_del(1)
#print "Kbc\n",K_i_H




