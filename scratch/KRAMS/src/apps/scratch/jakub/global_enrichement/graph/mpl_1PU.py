'''
Created on Apr 20, 2010

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

plt.axhline(1.)

#bilinear - 1PU
x1 = np.array([64,44,32,22,16,8,4], dtype = int)*2 
y1 = np.array([0.973495364189,0.96183681488,0.948169231415,0.926151096821,0.892445206642,0.816784799099,0.67709505558], dtype = float)
p1 = plt.semilogx(x1, y1, 'bo-') 

#biquadratic -lagrange 1PU
x2 = np.array([63,45,33,23,15,11,5], dtype = int)*2 
y2 = np.array([0.947873592377,0.927179100513,0.90982311964,0.868773639202,0.812013149261,0.761070668697,0.676301777363], dtype = float)
p2 = plt.semilogx(x2, y2, 'rD--') 

#bicubic -lagrange -1 PU
x3 = np.array([64,46,34,22,16,10,4], dtype = int)*2 
y3 = np.array([0.978135883808,0.968210220337,0.955681741238,0.927601695061,0.900899589062,0.81528031826,0.621305942535], dtype = float)
p3 = plt.semilogx(x3, y3, 'g^:') 

plt.xlim(7,150)
plt.ylim(0.60,1.1)
plt.legend( (p1[0], p2[0],p3[0]), (r'bilinear', r'biquadratic l.', r'bicubic l.'), loc='lower right', prop = {'size':20,
                                                                                  'family':'serif'})

#plt.xlabel(r'$\log(dof_x)$',fontsize = 20, family='serif')
#plt.ylabel(r'$\sigma_I^{max}$ [MPa]',fontsize = 20)


plt.annotate('exact solution', xy= (52,1.), xytext=(20,1.05),fontsize=18,family='serif',
            arrowprops=dict(facecolor='black',shrink=0.05)
            )


#plt.title('Plasticity: Return Mapping')


if __name__ == '__main__':
    plt.show()