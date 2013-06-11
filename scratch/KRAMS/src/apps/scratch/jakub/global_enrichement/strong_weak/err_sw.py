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

#plt.axhline(1.)

#bilinear
x1 = np.array([64,44,32,22,16,8,4], dtype = int)*2 
y1 = 1.-np.array([0.973495364189,0.96183681488,0.948169231415,0.926151096821,0.900899589062,0.816784799099,0.67709505558], dtype = float)
p1 = plt.loglog(x1, y1, 'bo-') 
#p1 = plt.loglog(x1, y1, 'bo-') 

#bilinear + weak
x2 = np.array([64,44,32,22,16,8,4], dtype = int)*2 
y2 = 1.-np.array([0.973972022533,0.962834477425,0.950031161308,0.930005013943,0.907991945744,0.842140555382, 0.74770539999], dtype = float)
p2 = plt.loglog(x2, y2, 'rD--') 
#p2 = plt.loglog(x2, y2, 'rD--') 

##bilinear + weak
#x3 = np.array([64,44,32,22,16,11,7], dtype = int)*2 
#y3 = 1.-np.array([0.973972022533,0.962834477425,0.950031161308,0.930005013943,0.907991945744,0.801617443562, 0.748394012451], dtype = float)
#p3 = plt.loglog(x3, y3, 'g^:') 


plt.xlim(7,150)
#plt.ylim(0.60,1.1)
plt.legend( (p1[0], p2[0]), (r'strong', r'strong+weak'), loc='upper right', prop = {'size':20})#,
                                                                                  #'family':'serif'})

plt.xlabel(r'DOF$_x$',fontsize = 20)#, family='serif')
plt.ylabel(r'maximum norm error',fontsize = 20)


#plt.annotate('exact solution', xy= (52,1.), xytext=(20,1.05),fontsize=18,family='serif',
#            arrowprops=dict(facecolor='black',shrink=0.05)
#            )


#plt.title('Plasticity: Return Mapping')


if __name__ == '__main__':
    plt.show()