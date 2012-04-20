'''
Created on 03.04.2012

@author: Axel
'''
from matplotlib import pyplot as plt
from numpy import log as ln, linspace
from math import pi as Pi
import numpy as np


def Ar(z, lf):
    #z is the dist between the 
    lf = np.float(lf)
    #######################
            #int#
    #######################
    
    if type(z) == int:
        if z == 0:
            return 1
        if z == lf / 2:
            return 0
        return lf * (lf - 2. * z + 2. * z * (.6931471806 + ln((1. / lf) * z))) / (lf - 2. * z) ** 2
    
    #######################
            #ARRAY#
    #######################
    
    res = 2 * z / lf * (ln(2 * z / lf) - 1) + 1
    if any(z == lf / 2):
        z = list(z)
        ind_h = z.index(lf / 2)
        res[ind_h] = 0
        z = np.array(z)
    if any(z == 0):
        z = list(z)
        ind_z = z.index(0) 
        res[ind_z] = 1 
        z = np.array(z)
    return res 

if __name__ == '__main__':
    z = linspace(0, 5, 400)
    plt.plot(z, Ar(z, 10))
    print Ar(z, 10)
    plt.show()
