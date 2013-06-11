'''
Created on Nov 24, 2009

@author: jakub
both phases clamped, traction prescribed 1. (/2 both phases), crack in matrix 90 deg
'''
import matplotlib
import numpy as np
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid.axislines import SubplotZero

matplotlib.rcParams['xtick.direction'] = 'out'
matplotlib.rcParams['ytick.direction'] = 'out'


# You can force all the contours to be the same color.
fig = plt.figure()
#
#plt.arrow(-400,0,800,0,head_width = 10)
#plt.arrow(0,-400,0,800,head_width = 10)
#plt.figtext(0.8,0.48,'$\sigma_1$',size=20,
#            horizontalalignment='center',
#            verticalalignment='top')
#plt.figtext(0.47,0.8,'$\sigma_2$',size=20,
#            horizontalalignment='left',
#            verticalalignment='center')

plt.axhline(1.)

#bilinear
x1 = np.array([64,44,32,22,16,8,4], dtype = int)*2 
y1 = np.array([0.923105638359,0.910130219959,0.894839259558,0.87002539402,0.841306947549,0.744013051223,0.571480475265], dtype = float)
p1 = plt.semilogx(x1, y1, 'bo-') 

#bicubic
x2 = np.array([64,46,34,22,16,10,4], dtype = int)*2 
y2 = np.array([0.928682253679,0.918037692774,0.902099521689,0.870117285041,0.838330181753,0.755096104431, 0.523177858396], dtype = float)
p2 = plt.semilogx(x2, y2, 'rD--') 

#bicubic
x3 = np.array([64,46,34,22,16,10,4], dtype = int)*2 
y3 = np.array([0.91847900712,0.905052613105,0.888192853473,0.852789402249,0.815622228873,0.737265600384, 0.503108964105], dtype = float)
p3 = plt.semilogx(x3, y3, 'g^:') 

plt.xlim(7,150)
plt.ylim(0.40,1.1)
plt.legend( (p1[0], p2[0],p3[0]), (r'bilinear', r'bicubic lagrange', r'bicubic serendipity'), loc='lower right', prop = {'size':20,
                                                                                  'family':'serif'})

#plt.xlabel(r'$\log(dof_x)$',fontsize = 20, family='serif')
#plt.ylabel(r'$\sigma_I^{max}$ [MPa]',fontsize = 20)


plt.annotate('exact solution', xy= (52,1.), xytext=(20,1.05),fontsize=18,family='serif',
            arrowprops=dict(facecolor='black',shrink=0.05)
            )


#plt.title('Plasticity: Return Mapping')


if __name__ == '__main__':
    plt.show()