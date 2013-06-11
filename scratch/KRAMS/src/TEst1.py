'''
Created on 14.12.2012

@author: acki
'''

from math import e
from matplotlib import pyplot as plt
import numpy as np
from math import pi as Pi
from scipy.stats import weibull_min
import scipy.stats as ss
from scipy.optimize import brentq

L_0 = 10.     #mm
eps_0 = 0.02   #
m = 4.          #Weibull modulus
tau = .1         #N/mm^2
r = 0.003        #mm radius fiber
sigma_app = 3000.   #N/mm^2
E_f = 100e3
eps_max = sigma_app / E_f
T = 2. * tau / r / E_f
dp = 100.          #Data Points
MCdp = 1000.
g_dp = 100.

eps_arr = np.linspace(1 / 1e6, 0.04, 1000)


def a(eps):
    return eps / T

def Pf_real(eps_f0):
    L1 = a(eps_f0) / 4
    return 1 - e ** (2 * (-a(eps_f0) * (eps_f0 / eps_0) ** m / (L_0 * (m + 1))) * (1 - (1 - L1 / a(eps_f0)) ** (m + 1)))
    #return 1 - e ** (-eps_f0 ** (m + 1) * E_f / eps_0 ** m * r / tau / L_0 / (m + 1) * (1 - (1 - L1 / a(eps_f0) ** (m + 1))))

def Pf_one_sided(eps_f0):
    L1 = a(eps_f0) / 4
    #return 1 - e ** (-(eps_f0 * (eps_f0 / eps_0) ** m + (-(L1 * T - eps_f0) / eps_0) ** m * L1 * T - (-(L1 * T - eps_f0) / eps_0) ** m * eps_f0) / T / (m + 1) / L_0)
    return 1 - e ** (2 * ((L1 * T - eps_f0) ** (m + 1) - eps_f0 ** (m + 1)) / (T * (m + 1) * L_0 * eps_0 ** m))
#1-exp(-(ef0*(ef0/eps0)^m+(-(L1*T-ef0)/eps0)^m*L1*T-(-(L1*T-ef0)/eps0)^m*ef0)

plt.plot(eps_arr, Pf_real(eps_arr), 'r')
plt.plot(eps_arr, Pf_one_sided(eps_arr), 'b')
plt.show()
