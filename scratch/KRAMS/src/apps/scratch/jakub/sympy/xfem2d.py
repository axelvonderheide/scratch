'''
Created on Sep 18, 2009

@author: jakub
'''
from sympy import symbols, Matrix, Plot, diff, sign,sin,abs, integrate,pprint

#from sympy.functions import sign

from sympy.functions.special.delta_functions import Heaviside, DiracDelta


Xi_jump = 0.5
x,y,x_i,Xi = symbols('x,y,x_i,Xi') 

def SGN_S(x):
    return x/abs(x)

def SGN_H(x):
    return 2 * Heaviside(x) - 1.

def Jump_S(x,x_i,Xi) :
    return 1./2.*(SGN_S(x-Xi)-SGN_S(x_i-Xi))

def Jump_H(x,x_i,Xi):
    return 1./2.*(SGN_H(x-Xi)-SGN_H(x_i-Xi))

N_mtx_S = Matrix(1,4,[(x-1)*(y-1),-x*(y-1), (x-1)*(y-1)*Jump_S(x,0.,Xi_jump), -x*(y-1)*Jump_S(x,1.,Xi_jump)])
N_mtx_H = Matrix(1,4,[(x-1)*(y-1),-x*(y-1), (x-1)*(y-1)*Jump_H(x,0.,Xi_jump), -x*(y-1)*Jump_H(x,1.,Xi_jump)])

print "N_s\n",N_mtx_S
print "N_h\n",N_mtx_H

dx = lambda a: diff(a,x)
dy = lambda a: diff(a,y)

#p = Plot()
#for i,N in enumerate(N_mtx_S):
#    p[i]=N,[x,0,1,100]
    
B_mtx_S_x = N_mtx_S.applyfunc(dx)
B_mtx_H_x = N_mtx_H.applyfunc(dx)

#print 'B_sx\n',B_mtx_S_x
#print 'B_hx\n',B_mtx_H_x

B_mtx_S_xy = N_mtx_S.applyfunc(dy)
B_mtx_H_xy = N_mtx_H.applyfunc(dy)

#print 'B_sxy\n',B_mtx_S_xy
#print 'B_hxy\n',B_mtx_H_xy

B_mtx_S = B_mtx_S_x.col_join(B_mtx_S_xy)
B_mtx_H = B_mtx_H_x.col_join(B_mtx_H_xy)

print 'B_s\n',B_mtx_S
print 'B_h\n',B_mtx_H

B_F = Matrix(2,4,[-1 + y, 1 - y, 0, 0,\
                  -1 + x, -x, 0, 0])
                  
B_H = Matrix(2,4,[0, 0, -(1 - y)*Heaviside(-0.5 + x) , -(1 - y)*(1.0 - Heaviside(-0.5 + x)),\
                  0, 0, -(1 - x)*Heaviside(-0.5 + x), x*(1.0 - Heaviside(-0.5 + x))])

B_D = Matrix(2,4,[0, 0, (1 - x)*(1 - y)*DiracDelta(-0.5 + x), x*(1 - y)*DiracDelta(-0.5 + x),\
                  0, 0, 0, 0])

Ef, Em, Es = symbols('Ef,Em,Es') 
Es = Ef #+ Em *(1 - DiracDelta(Xi_jump))
D_mtx = Matrix([[Es,0],\
                [0,(Es)/2.]])#nu=0.

#K_r_S = B_mtx_S.T * D_mtx*B_mtx_S
#K_r_H = B_mtx_H.T * D_mtx*B_mtx_H
#
#K_r_Sdata = K_r_S.subs(Ef,20).subs(Em,20)
#K_r_Hdata = K_r_H.subs(Ef,20).subs(Em,20)

#print 'K_r_s\n',K_r_S.shape
#print 'K_r_h\n',K_r_H.shape

integrx = lambda a:integrate(a,(x,0,1))
integrx_l = lambda a:integrate(a,(x,0.,0.5))
integrx_r = lambda a:integrate(a,(x,0.5,1.))

integry = lambda a:integrate(a,(y,0,1))

#K_i_H = K_r_Hdata.applyfunc(integry).applyfunc(integrx)
#print "KiH\n",K_i_H
#
#K_i_Hc = K_r_Hdata.applyfunc(integry).applyfunc(integrx_l)\
#        +K_r_Hdata.applyfunc(integry).applyfunc(integrx_r)
#print "KiHc\n",K_i_Hc

#print "diff\n",K_i_H-K_i_Hc

K_F = ((B_F.T * D_mtx*B_F).subs(Ef,20).subs(Em,20)).applyfunc(integry).applyfunc(integrx)
print "KF\n",K_F

K_H = ((B_H.T * D_mtx*B_H).subs(Ef,20).subs(Em,20)).applyfunc(integry).applyfunc(integrx_l)
print "KH\n",K_H