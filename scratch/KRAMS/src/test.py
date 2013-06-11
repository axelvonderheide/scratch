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
MCdp = 1.
g_dp = 50.

def H(x):
    return x >= 0

def a(eps):
    return r * eps * E_f / (2 * tau)

def res_eps(z, eps_app, eps_strength):
    #plt.plot(z, eps_app - T * z)
    #plt.plot(z, np.ones(len(z)) * eps_strength)
    #plt.plot(z, eps_strength - (eps_app - T * z))
    #plt.show()
    return eps_strength - (eps_app - T * z)

def break_gen(z, eps_app):
        #random numbers
        r1 = np.random.rand(dp)
        r2 = np.random.rand(dp)
        
        #picking lowest numbers
        r_low = np.maximum(r1, r2)

        #PPF of epsilon fed with random numbers between 0-1
        delta_z = a(eps_app) / (dp - 1)
        eps_strength = eps_0 * (-np.log(r_low) * L_0 / delta_z) ** (1 / m)
        
        
        #####Mirroring Distr
        x1_MDISTR = np.linspace(0.001, 1, 20)
        x2_MDISTR = np.linspace(0.000001, 0.05, 20)
        y1 = eps_0 * (-np.log(1 - x1_MDISTR) * L_0 / delta_z) ** (1 / m)
        #print (np.log(1 - x1_MDISTR) * L_0 / delta_z)
        y2 = 1 - e ** (-delta_z / L_0 * (x2_MDISTR / eps_0) ** m)
        
        plt.plot(x2_MDISTR, y2, 'b')
        plt.plot(y1, x1_MDISTR, 'k')
        plt.show()
        
        
        #########################
        
        
        
        #test of strength and strain
        broken_test = res_eps(z, eps_app, eps_strength)
        
        #picking lowest number
        min_v = np.min(broken_test)
        return min_v < 0


        
        
def P_f_func(eps):
            t = -2. *a(eps) / L_0 * (eps / eps_0) ** m / (m + 1.)
            return 1 - e ** (t)


   
mean_gen_list = []  
def P_f(eps):
    z = np.linspace(0, a(eps), dp)
    
    eps_MC = np.linspace(eps, eps, MCdp)
    for eps_i in eps_MC:
        nob = break_gen(z, eps_i)
        mean_gen_list.append(nob)
    return np.mean(mean_gen_list)


P_f_list = []
def MC_Test():
    
    eps_arr = np.linspace(eps_max / 10000, eps_max, g_dp)
    
    for eps in eps_arr:
        P_f_list.append(P_f(eps))

    plt.plot(eps_arr, P_f_list)
    plt.plot(eps_arr, P_f_func(eps_arr))
    plt.show()
    
    
MC_Test()




