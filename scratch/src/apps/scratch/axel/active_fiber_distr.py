'''
Created on 28.03.2012

@author: Axel
'''
from stats.pdistrib.sinus_distribution import sin_distr
from math import pi as Pi
import numpy as np

l = 10

phi = np.linspace(0, Pi / 2, 50)
x = np.linspace(-l / 2, l / 2, 50)
y = np.linspace(-l / 2, l / 2, 50)

def H(x):
    return np.sign(np.sign(x) + 1.)

x_pdf = l / 2 #Def: [-l/2;l/2]
Ic = H((1 / 2) * l * np.cos(phi) - abs(x) - abs(y))



