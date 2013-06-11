#from math import e
from matplotlib import pyplot as plt
import numpy as np
from math import pi as Pi
from scipy.stats import weibull_min
import scipy.stats as ss
from scipy.optimize import brentq
from scipy.interpolate import interp1d


#Data Points
dp = 250.         
MCdp = 50000.
g_dp = 5.
dp_ana = 300.


#Parameters

L_0 = 100.     #mm
eps_0 = 0.02   #
s = eps_0
m = 4.          #Weibull modulus
tau = .1         #N/mm^2
r = 0.003        #mm radius fiber
sigma_app = 4000.   #N/mm^2
E_f = 100e3
eps_max = sigma_app / E_f
T = 2. * tau / (r * E_f)


#Code

def H(x):
    return x >= 0

def a(eps):
    return eps / T


#Monte Carlo

def res(eps_app, z, eps_strength):
    return np.min(eps_strength - (eps_app * np.ones(len(z)) - T * z))

def res_eps(eps_app, z, eps_strength):
    #plt.plot(z, eps_app - T * z)
    #plt.plot(z, np.ones(len(z)) * eps_strength)
    #plt.plot(z, eps_strength - (eps_app - T * z))
    #plt.show()
    return eps_strength - (eps_app * np.ones(len(z)) - T * z)


def z_gen(z, eps_app, eps_strength):
        #find out eps, that prevails at the break of the fiber
        try:
            eps_br = brentq(res, 0., eps_app, args=(z, eps_strength))
            z_index = np.argmin(res_eps(eps_br, z, eps_strength))
            z_pullout = z[z_index]
            return z_pullout, 1
    
        except:
            return 0, 0
            
        
        

 
def z_mean(eps, eps_strength):
    
    mean_gen_list = [] 
    
    z = np.linspace(0, a(eps), dp)
    
    eps_MC = np.linspace(eps, eps, MCdp)
    
    for eps_i in eps_MC:
        z_i, bool = z_gen(z, eps_i, eps_strength)
        if bool == 1:
            mean_gen_list.append(z_i)
    
    
    if len(mean_gen_list) == 0:
        mean_gen_list.append(0)   
        
           
    return np.mean(mean_gen_list)

def gen_eps_strength(z_main):
        #random numbers
        r1 = np.random.rand(g_dp)
        r2 = np.random.rand(g_dp)
        #picking lowest numbers
        r_low = np.maximum(r1, r2)
        #PPF of epsilon fed with random numbers between 0-1
        delta_z = z_main[-1] / (len(z_main) - 1)
        return eps_0 * (-np.log(r_low) * L_0 / delta_z) ** (1 / m)

########################################


############# analytical #################

def P_f(eps, L1):
            t = (2 * -a(eps) * (eps / eps_0) ** m / (L_0 * (m + 1)) * (1 - (1 - L1 / a(eps)) ** (m + 1)))
            return 1. - np.exp(t)
        
def z_mean_ana(eps):
            z = np.linspace(0, a(eps), dp_ana)
            a_eps = a(eps)
            int_z = np.trapz(P_f(eps, z), z)
            return a_eps - int_z / P_f(eps, a_eps)

#########################################

def CDFa(e):
    T = 2. * tau / r / E_f
    a = e / T
    return 1. - np.exp(-a * 2.* (e / s) ** m / (m + 1) / L_0)

def PDFa(e):
    T = 2. * tau / r / E_f
    a = e / T
    return np.exp(-a * 2.* (e / s) ** m / L_0 / (m + 1)) * 2.*(e / s) ** m / (T * L_0)

def CDFL(e, L):
    T = 2. * tau / r / E_f
    a = e / T
    return 1. - np.exp(-a * 2.* (e / s) ** m * (1 - (1 - L / a) ** (m + 1)) / (m + 1) / L_0)

def scalar_mu_L(e):
    c = 1. / CDFa(e)
    T = 2. * tau / r / E_f
    a = e / T
    z_arr = np.linspace(0, a, dp_ana)
    integ = np.trapz(CDFL(e, z_arr), z_arr)
    z_arr * PDFa(e)
    return a - c * integ

def mu_L(e_arr):
    muL = []
    for e in e_arr:
        muL.append(scalar_mu_L(e))
    return np.array(muL)

def scalar_mu_L_rostar(e):
    e_arr = np.linspace(0.0001, e, dp_ana)
    muL_arr = mu_L(e_arr)
    PDF = PDFa(e_arr)
    c = 1. / CDFa(e)
    integ = np.trapz(PDF * muL_arr, e_arr)
    return c * integ

def mu_L_rostar(e_arr):
    muL = []
    for e in e_arr:
        print np.round(e / (e_arr[-1] - e_arr[0]) * 100, 2), '% Loading  analytic curve'
        muL.append(scalar_mu_L_rostar(e))
    return np.array(muL)


#############################################################

z_ana_list = []
z_mean_list = []
def MC_Test():
    
    eps_strength = gen_eps_strength()
    #eps_arr = np.linspace(0.02, eps_max, g_dp)
    eps_arr = np.linspace(1e-6, 0.035, g_dp)
    z_main = eps_arr / T
    for eps in eps_arr:
        print np.round(eps / (eps_arr[-1] - eps_arr[0]) * 100, 2), '% LOADING MC Simulation'  
        z_mean_list.append(z_mean(eps))
        z_ana_list.append(z_mean_ana(eps))
        
    plt.plot(eps_arr, mu_L_rostar(eps_arr), label='mu_L Rypl')
    ###########
    print 'DONE'
    ###########
    plt.xlabel('eps crack')
    plt.ylabel('mu_L')
    plt.legend(loc='best')
    z_ana_list[0] = 0
    plt.plot(eps_arr, z_mean_list, 'r', label='z_mean')
    plt.plot(eps_arr, z_ana_list, 'k')
    plt.show()
    
    
MC_Test()
