'''
Created on Apr 27, 2010

@author: jakub
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

#plt.axhline(1.)
e_exact = 0.58476470554775783

#bilinear - 1PU - strong
#x1 = np.array([64,44,32,22,16,8,4], dtype = int)*2 
#y1 = np.array([0.02355117,0.03457074,0.04805082,0.07113152,0.09981872,0.213483019308,0.465710846741], dtype = float)
#two values dropped
x1 = np.array([64,32,16,8,4], dtype = int)*2 
y1 = np.array([ 0.03438816, 0.06776366, 0.13120908,0.2430208, 0.39533837], dtype = float)/e_exact
p1 = plt.loglog(x1, y1, 'bo-') 

#bilinear - 1PU - strong+weak
x2 = np.array([64,32,16,8,4], dtype = int)*2 
y2 = np.array([ 0.03438816, 0.06776366, 0.13120905,0.24301965, 0.39527241], dtype = float)/e_exact
p2 = plt.loglog(x2, y2, 'rD--') 

plt.xlim(8,150)
plt.ylim(1.e-3,9.e-1)
#plt.legend( (p1[0], p2[0],p3[0]), (r'bilinear', r'biquadratic l.', r'bicubic l.'), loc='upper right', prop = {'size':20,
#                                                                                  'family':'serif'})

#plt.xlabel(r'$\log(dof_x)$',fontsize = 20, family='serif')
#plt.ylabel(r'$\sigma_I^{max}$ [MPa]',fontsize = 20)


#plt.annotate('exact solution', xy= (52,1.), xytext=(20,1.05),fontsize=18,family='serif',
#            arrowprops=dict(facecolor='black',shrink=0.05)
#            )


#plt.title('Plasticity: Return Mapping')


if __name__ == '__main__':
    plt.show()