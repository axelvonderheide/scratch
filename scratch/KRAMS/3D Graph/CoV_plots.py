'''
Created on 2.2.2013

@author: Q
'''
from indep_CB_model import CBResidual
from stats.spirrid.spirrid import SPIRRID
from stats.spirrid.rv import RV
from matplotlib import pyplot as plt
import numpy as np
from scipy.optimize import fmin
from etsproxy.mayavi import mlab
from stats.spirrid import make_ogrid as orthogonalize
from scipy.optimize import fsolve
from math import pi
from scipy.special import gammainc, gamma
def H(x):
    return x >= 0.0

def get_scale(mu_xi, m, tau, r):
        def optimize(s):
            p = np.linspace(0., .9999, 1000)
            T = 2. * tau / r
            ef0_break = (-0.5 * np.log(1. - p) * T / E_f * (m + 1) * s ** m) ** (1. / (m + 1))
            return np.trapz(1 - p, ef0_break) - mu_xi
        return fsolve(optimize, mu_xi)
    
def rand_tau_r(E_f, V_f, mu_tau, mu_r, mu_xi, m_list, CoV_tau_range, CoV_r_range, dp_tau, dp_r, sV0):
    #rand tau and rand r
    for mi in m_list:
        s = get_scale(mu_xi, mi, mu_tau, mu_r)
        sV0 = float(s * (pi * mu_r ** 2) ** (1. / mi))
        Pf = RV('uniform', loc=0.0, scale=1.0)
        w_arr = np.linspace(0, 1.2, 30)
        cb = CBResidual(include_pullout=True)
        # loc scale generation for specific CoV
        #CoVtau
        CoV_tau_arr = np.linspace(CoV_tau_range[0], CoV_tau_range[1], dp_tau)
        loc_tau_arr = mu_tau - CoV_tau_arr * mu_tau * 3 ** 0.5
        scale_tau_arr = 2 * mu_tau - 2 * loc_tau_arr
        #CoVr
        CoV_r_arr = np.linspace(CoV_r_range[0], CoV_r_range[1], dp_r)
        loc_r_arr = mu_r - CoV_r_arr * mu_r * 3 ** 0.5
        scale_r_arr = 2 * mu_r - 2 * loc_r_arr
        #shaping for mayavi
        e_arr = orthogonalize([CoV_tau_arr, CoV_r_arr])
        x_axis = e_arr[0]
        y_axis = e_arr[1]
        #TAU gen Tuple of [loc,scale]
        stats_tau = []
        for s in range(dp_tau):
            stats_tau.append(RV('uniform', loc=loc_tau_arr[s], scale=scale_tau_arr[s]))
    
        stats_r = []
        for s in range(dp_r):
            stats_r.append(RV('uniform', loc=loc_r_arr[s], scale=scale_r_arr[s]))
             
        #r gen Tuple of [loc,scale]
        
        res_array = np.zeros((dp_tau, dp_r))
        for i_tau, taui in enumerate(stats_tau): 
            
                for i_r, ri in enumerate(stats_r):
                    
                    total = SPIRRID(q=cb,
                            sampling_type='MCS',
                            evars=dict(w=w_arr),
                            tvars=dict(tau=taui, E_f=E_f, V_f=V_f, r=ri,
                                       m=mi, sV0=sV0, Pf=Pf),
                            n_int=70)
                    if isinstance(ri, RV):
                        r_arr = np.linspace(ri.ppf(0.001), ri.ppf(0.999), 200)
                        Er = np.trapz(r_arr ** 2 * ri.pdf(r_arr), r_arr)
                    else:
                        Er = ri ** 2
                    result = total.mu_q_arr / Er
                    
                    sigma_c = np.max(result)
                    
                    if sigma_c == result[-1]:
                        print "w_arr too short"
                        
                    res_array[i_tau, i_r] = sigma_c
        #mayaviplot
        
        mlab.surf(x_axis, y_axis, res_array)
        #mlab.xlabel("rand tau")
        #mlab.ylabel("rand r")
        #mlab.zlabel("sigma")
    mlab.show()
    
###########################################################################################################

###########################################################################################################
"""
def det_tau_r(E_f, V_f, m_list, mu_xi, dp_tau, dp_r, tau_range, r_range):
    #loop with rand tau and rand r
    for mi in m_list:
        Pf = RV('uniform', loc=0.0, scale=1.0)
        w_arr = np.linspace(0, 40, 300)
        cb = CBResidual(include_pullout=True)
        tau_arr = np.linspace(tau_range[0], tau_range[1], dp_tau)
        r_arr = np.linspace(r_range[0], r_range[1], dp_r)
        e_arr = orthogonalize([tau_arr, r_arr])
        x_axis = e_arr[0]
        y_axis = e_arr[1]
        res_array = np.zeros((dp_tau, dp_r))
        for i_tau, taui in enumerate(tau_arr):
                for i_r, ri in enumerate(r_arr):
                    s = get_scale(mu_xi, mi, taui, ri)
                    sV0 = float(s * (pi * ri ** 2) ** (1. / mi))
                    total = SPIRRID(q=cb,
                            sampling_type='MCS',
                            evars=dict(w=w_arr),
                            tvars=dict(tau=taui, E_f=E_f, V_f=V_f, r=ri,
                                       m=mi, sV0=sV0, Pf=Pf),
                            n_int=10000)
                    if isinstance(ri, RV):
                        r_arr = np.linspace(ri.ppf(0.001), ri.ppf(0.999), 200)
                        Er = np.trapz(r_arr ** 2 * ri.pdf(r_arr), r_arr)
                    else:
                        Er = ri ** 2
                    result = total.mu_q_arr / Er
                    sigma_c = np.max(result)
                    #plt.plot(w_arr, result)
                    #plt.show()
                    if sigma_c == result[-1]:
                        print "w_arr too short"
                        pass
                    res_array[i_tau, i_r] = sigma_c
    
        #mayaviplot
        mlab.surf(x_axis, 10 * y_axis , res_array)
        #mlab.xlabel("det tau")
        #mlab.ylabel("det r")
        #mlab.zlabel("sigma")
 
    tau_ax = np.array([0.0000001, tau_range[1] * 1.2])
    r_ax = np.array([0.0000001, r_range[1] * 12])
    z_ax = np.zeros((2, 2))
    z_ax[0, 0] = 20.
    e_arr = orthogonalize([tau_ax, r_ax])
    x_ax = e_arr[0]
    y_ax = e_arr[1]
    mlab.surf(x_ax, y_ax, z_ax, representation='points')
    mlab.view(0., 0.)
    mlab.xlabel("det tau")
    mlab.ylabel("det r")
    mlab.zlabel("sigma")
    mlab.show()
""" 

def det_tau_r(E_f, V_f, m_list, mu_xi, dp_tau, dp_r, tau_range, r_range, sV0):
    w_arr = np.linspace(0, 20, 300)
    tau_arr = np.linspace(tau_range[0], tau_range[1], dp_tau)
    r_arr = np.linspace(r_range[0], r_range[1], dp_tau)
    e_arr = orthogonalize([tau_arr, r_arr])
    x_axis = e_arr[0]
    y_axis = e_arr[1]
    for mi in m_list:
            res_array = np.zeros((dp_tau, dp_r))
            for i_tau, taui in enumerate(tau_arr):
                    for i_r, ri in enumerate(r_arr):
                        s = get_scale(mu_xi, mi, taui, ri)
                        print s
                        #sV0 = float(s * (pi * mu_r ** 2) ** (1. / mi))
                        T = 2. * taui / ri
                        s0 = ((T * (mi + 1) * sV0 ** mi) / (2. * E_f * pi * ri ** 2)) ** (1. / (mi + 1))
                        print s0
                        k = np.sqrt(T / E_f)
                        ef0 = k * np.sqrt(w_arr)
                        G = 1 - np.exp(-(ef0 / s0) ** (mi + 1))
                        mu_int = ef0 * E_f * V_f * (1 - G)
                        I = s0 * gamma(1 + 1. / (mi + 1)) * gammainc(1 + 1. / (mi + 1), (ef0 / s0) ** (mi + 1))
                        mu_broken = E_f * V_f * I / (mi + 1)
                        result = mu_int + mu_broken
                        sigma_c = np.max(result)
                        if sigma_c == result[-1]:
                            print "w_arr too short"
                            pass
                        res_array[i_tau, i_r] = sigma_c
                    
            mlab.surf(x_axis, y_axis, res_array / 100.)
    mlab.xlabel("det tau")
    mlab.ylabel("det r")
    mlab.zlabel("sigma")
    mlab.show()
    
#########################################################################

######################################################################## 

def rand_tau_det_r(E_f, V_f, CoV_tau_range, r_range, mu_tau , dp_tau, dp_r, sV0):
    for mi in m_list:
        #loop with rand tau and rand r
        Pf = RV('uniform', loc=0.0, scale=1.0)
        w_arr = np.linspace(0, 5, 300)
        cb = CBResidual(include_pullout=True)
        # CoV generation
        CoV_tau_arr = np.linspace(CoV_tau_range[0], CoV_tau_range[1], dp_tau)
        loc_tau_arr = mu_tau - CoV_tau_arr * mu_tau * 3 ** 0.5
        scale_tau_arr = 2 * mu_tau - 2 * loc_tau_arr
        r_arr = np.linspace(r_range[0], r_range[1], dp_r)
        e_arr = orthogonalize([CoV_tau_arr, r_arr])
        x_axis = e_arr[0]
        y_axis = e_arr[1]
        #TAU gen Tuple of [loc,scale]
        stats_tau = []
        for s in range(dp_tau):
            stats_tau.append(RV('uniform', loc=loc_tau_arr[s], scale=scale_tau_arr[s]))
    
        #r gen Tuple of [loc,scale]
        res_array = np.zeros((dp_tau, dp_r))
        for i_tau, taui in enumerate(stats_tau): 
                for i_r, ri in enumerate(r_arr):
                    s = get_scale(mu_xi, mi, mu_tau, ri)
                    sV0 = float(s * (pi * ri ** 2) ** (1. / mi))
                    total = SPIRRID(q=cb,
                            sampling_type='MCS',
                            evars=dict(w=w_arr),
                            tvars=dict(tau=taui, E_f=E_f, V_f=V_f, r=ri,
                                       m=mi, sV0=sV0, Pf=Pf),
                            n_int=60)
                    if isinstance(ri, RV):
                        r_arr = np.linspace(ri.ppf(0.001), ri.ppf(0.999), 200)
                        Er = np.trapz(r_arr ** 2 * ri.pdf(r_arr), r_arr)
                    else:
                        Er = ri ** 2
                    result = total.mu_q_arr / Er
                    sigma_c = np.max(result)
                    if sigma_c == result[-1]:
                        print "w_arr too short"
                    res_array[i_tau, i_r] = sigma_c
        #mayaviplot
        mlab.surf(x_axis, y_axis * 50, res_array, warpscale=0.1)
        mlab.view(0., 0.)
        mlab.xlabel("rand tau")
        mlab.ylabel("det r")
        mlab.zlabel("sigma")
    mlab.show()

##############################################
def rand_r_det_tau(E_f, V_f, tau_range, CoV_r_range, mu_r , dp_tau, dp_r, sV0):
    #loop with rand tau and rand r
    for mi in m_list:
        Pf = RV('uniform', loc=0.0, scale=1.0)
        w_arr = np.linspace(0, 2.0, 30)
        cb = CBResidual(include_pullout=True)
        tau_arr = np.linspace(tau_range[0], tau_range[1], dp_tau)
    
        CoV_r_arr = np.linspace(0.0001, 0.5, dp_r)
        loc_r_arr = mu_r - CoV_r_arr * mu_r * 3 ** 0.5
        scale_r_arr = 2 * mu_r - 2 * loc_r_arr
        e_arr = orthogonalize([tau_arr, CoV_r_arr])
        x_axis = e_arr[0]
        y_axis = e_arr[1]
        stats_r = []
        for s in range(dp_r):
            stats_r.append(RV('uniform', loc=loc_r_arr[s], scale=scale_r_arr[s]))
    
        #gen Tuple of [loc,scale]
    
        res_array = np.zeros((dp_tau, dp_r))
        for i_tau, taui in enumerate(tau_arr):
    
                for i_r, ri in enumerate(stats_r):
                    s = get_scale(mu_xi, mi, taui, mu_r)
                    sV0 = float(s * (pi * mu_r ** 2) ** (1. / mi))
                    total = SPIRRID(q=cb,
                            sampling_type='MCS',
                            evars=dict(w=w_arr),
                            tvars=dict(tau=taui, E_f=E_f, V_f=V_f, r=ri,
                                       m=mi, sV0=sV0, Pf=Pf),
                            n_int=60)
                    if isinstance(ri, RV):
                        r_arr = np.linspace(ri.ppf(0.001), ri.ppf(0.999), 200)
                        Er = np.trapz(r_arr ** 2 * ri.pdf(r_arr), r_arr)
                    else:
                        Er = ri ** 2
                    result = total.mu_q_arr / Er
                    sigma_c = np.max(result)
                    if sigma_c == result[-1]:
                        print "w_arr too short"
                        pass
                    res_array[i_tau, i_r] = sigma_c
    
        #mayaviplot
        mlab.surf(x_axis, y_axis, res_array)
        #
        mlab.view(0., 0.)
        mlab.xlabel("det tau")
        mlab.ylabel("rand r")
        mlab.zlabel("sigma")
    mlab.show()



#Parameter


#PLOTS
#########
"Params"
#########
E_f, V_f, sV0 = 200e3, 0.01, 3.e-3

#########################
"random tau and random r"
#########################

#mu of rndm variables
mu_tau = 0.1
mu_r = 0.01
mu_xi = .008

#plots of different m's with fix mu_xi
m_list = [5., 20.]

#Ranges of CoV's
CoV_tau_range = [0.00001, 0.5]
CoV_r_range = [0.00001, 0.5]
dp_tau = 5
dp_r = 5

#Plot
#rand_tau_r(E_f, V_f, mu_tau, mu_r, mu_xi, m_list, CoV_tau_range, CoV_r_range, dp_tau, dp_r,sV0)

###################
"det tau and det r"
###################

#ranges of tau and r
tau_range = [.05, 0.2]
r_range = [.001, .02]

#xi params
m_list = [ 5.]
mu_xi = .008
#Datapoints
dp_tau = 8
dp_r = 8
#Plot
det_tau_r(E_f, V_f, m_list, mu_xi, dp_tau, dp_r, tau_range, r_range, sV0)

###################
"rand tau and det r"
###################
#mu of rndm variables
#plots of different m's with fix mu_xi
mu_xi = .008
m_list = [5., 20.]

#Ranges of variables
mu_tau = 0.1
CoV_tau_range = [0.00001, 0.5]
r_range = [0.01, 0.03]
#Data Points
dp_tau = 5
dp_r = 10
#Plot
#rand_tau_det_r(E_f, V_f, CoV_tau_range, r_range, mu_tau , dp_tau, dp_r,sV0)

###################
"det tau and rand r"
###################
#mu of rndm variables
#plots of different m's with fix mu_xi
mu_xi = .008
m_list = [5., 20.]

#Ranges of variables
mu_r = 0.01
CoV_r_range = [0.00001, 0.5]
tau_range = [0.05, 0.15]
#Data Points
dp_tau = 5
dp_r = 10
#Plot
#rand_r_det_tau(E_f, V_f, tau_range, CoV_r_range, mu_r , dp_tau, dp_r,sV0)






