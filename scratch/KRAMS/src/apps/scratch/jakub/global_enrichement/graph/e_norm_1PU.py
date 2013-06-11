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

#bilinear - 1PU
#x1 = np.array([64,44,32,22,16,8,4], dtype = int)*2 
#y1 = np.array([0.02355117,0.03457074,0.04805082,0.07113152,0.09981872,0.213483019308,0.465710846741], dtype = float)
#two values dropped
x1 = np.array([64,32,16,8,4], dtype = int)*2 
y1 = np.array([ 0.03438816, 0.06776366, 0.13120908,0.2430208, 0.39533837], dtype = float)/e_exact
p1 = plt.loglog(x1, y1, 'bo-') 

#biquadratic -lagrange 1PU
x2 = np.array([63,46,33,23,15,11,5], dtype = int)*2 
y2 = np.array([0.0025266,0.00493783, 0.00815197,0.01877062, 0.04366996,0.08012097,0.19033243], dtype = float)/e_exact
p2 = plt.loglog(x2, y2, 'rD--') 

#biquadratic -lagrange 1PU
x2 = np.array([63,46,33,23,15,11,5], dtype = int)*2 
y2 = np.array([0.0025266,0.00493783, 0.00815197,0.01877062, 0.04366996,0.08012097,0.19033243], dtype = float)/e_exact
p2 = plt.loglog(x2, y2, 'rD--') 
##
#bicubic -lagrange -1 PU
x3 = np.array([64,34,22,16,10,4], dtype = int)*2 
y3 = np.array([ 0.00106053, 0.00389852,0.00971986,0.01916438, 0.05284005,0.41319661], dtype = float)/e_exact
p3 = plt.loglog(x3, y3, 'g^:') 

##bicubic -lagrange -3 PU
#x3 = np.array([64,34,22,16,10,4], dtype = int)*2 
#y3 = np.array([ 0.00106053, 0.00389852,0.00971986,0.01916438, 0.05284005,0.41319661], dtype = float)/e_exact
#p3 = plt.loglog(x3, y3, 'm+:') 

##bicubic -serendipity -1 PU
#x4 = np.array([64,34,22,16,10,4], dtype = int)*2 
#y4 = np.array([ 0.00098815, 0.00360036,0.0089472, 0.01759968, 0.04894311,0.41200072], dtype = float)/e_exact
#p4 = plt.loglog(x3, y3, 'ks-.') 



#biquartic -lagrange -1 PU
#x4 = np.array([53,45,37,29,21,13,5], dtype = int)*2 
#y4 = np.array([0.00571694,0.00662573, 0.00880001,0.01103427, 0.01730021,0.05850203,0.40842632], dtype = float)/e_exact
#p4 = plt.loglog(x4, y4, 'ks-.') 
#
plt.xlim(8,130)
plt.ylim(1.5e-3,9.e-1)
#plt.legend( (p1[0], p2[0],p3[0]), (r'bilinear', r'biquadratic l.', r'bicubic l.'), loc='upper right', prop = {'size':20,
#                                                                                  'family':'serif'})
plt.legend( (p1[0], p2[0],p3[0]), (r'Q4', r'Q9', r'Q16'), loc='upper right', prop = {'size':20})

plt.xlabel(r'DOF$_x$',fontsize = 20)#, family='serif')
plt.ylabel(r'energy norm error',fontsize = 20)


#plt.annotate('exact solution', xy= (52,1.), xytext=(20,1.05),fontsize=18,family='serif',
#            arrowprops=dict(facecolor='black',shrink=0.05)
#            )


#plt.title('Plasticity: Return Mapping')


if __name__ == '__main__':
    plt.show()