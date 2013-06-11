'''
Created on Apr 26, 2010

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

#bilinear - 1PU
#x1 = np.array([64,44,32,22,16,8,4], dtype = int)*2 
#y1 = np.array([0.02355117,0.03457074,0.04805082,0.07113152,0.09981872,0.213483019308,0.465710846741], dtype = float)
#two values dropped
x1 = np.array([64,32,16,8,4], dtype = int)*2 
y1 = np.array([ 6.49244685e-05, 0.00025622,0.0009884,0.00362505, 0.01159784], dtype = float)
p1 = plt.loglog(x1, y1, 'bo-') 

##biquadratic -lagrange 1PU
#x2 = np.array([63,33,15,11,5], dtype = int)*2 
#y2 = np.array([0.04727703,0.08567322,0.20382086,0.28105889,0.44963877], dtype = float)
#p2 = plt.loglog(x2, y2, 'rD--') 
##
##bicubic -lagrange -1 PU
#x3 = np.array([64,34,16,10], dtype = int)*2 
#y3 = np.array([0.978135883808,0.955681741238,0.900899589062,0.43784462], dtype = float)
#p3 = plt.loglog(x3, y3, 'g^:') 
#
plt.xlim(5,200)
plt.ylim(2.e-5,2.e-2)
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