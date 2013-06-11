from math import e
from matplotlib import pyplot as plt
import numpy as np
from math import pi as Pi
   


L_0 = 100.     #mm
eps_0 = 0.02   #
m = 4.          #Weibull modulus
tau = .1         #N/mm^2
r = 0.003        #mm radius fiber
sigma_app = 5000.   #N/mm^2
E_f = 100e3

eps_max = sigma_app / E_f
T = 2. * tau / r / E_f
dp = 100.          #Data Points


def a(eps):
    return r * eps * E_f / (2 * tau)


def P_f_cum(eps_f0):
        z = np.linspace(0, a(eps_f0), dp)
        exp_term = -2. *(a(eps_f0) / (dp - 1.)) / L_0 * ((eps_f0 - T * z) / eps_0) ** m 
        exp_term[0] = -1. *(a(eps_f0) / (dp - 1.)) / L_0 * ((eps_f0) / eps_0) ** m 
        exp_cum = np.cumsum(exp_term)
        return 1 - e ** exp_cum
        
def P_f_func(eps):
            t = -2. *a(eps) / L_0 * (eps / eps_0) ** m / (m + 1.)
            return 1 - e ** (t)
        
P_f_cum_list = []
def P_f_loop():      
    for eps_f0 in eps_arr:
        P_f = P_f_cum(eps_f0)[-1]
        P_f_cum_list.append(P_f)
    return 0


eps_arr = np.linspace(0, sigma_app / E_f, dp)
P_f_loop()
plt.plot(P_f_cum_list, eps_arr)
plt.plot(P_f_func(eps_arr), eps_arr)
#print 1 - P_f_func(eps_arr)[-1] / P_f_cum(z1)[-1]
#plt.plot(z1, eps(z1))
plt.show()
