'''
Created on 2.2.2013

@author: Q
'''
from indep_CB_model import CBResidual
from stats.spirrid.spirrid import SPIRRID
from stats.spirrid.rv import RV
from matplotlib import pyplot as plt
import numpy as np

def H(x):
    return x >= 0.0

def general_diagram():
    E_f, V_f, r, tau, m, sV0 = 200e3, 0.01, 0.01, .1, 7., 3.e-3
    Pf = RV('uniform', loc=0.0, scale=1.0)
    w_arr = np.linspace(0,1.2,200)
    cb = CBResidual(include_pullout=True)
    total = SPIRRID(q=cb,
                sampling_type='TGrid',
                evars=dict(w=w_arr),
                tvars=dict(tau=tau, E_f=E_f, V_f=V_f, r=r,
                           m=m, sV0=sV0, Pf=Pf),
                n_int=1000)
    if isinstance(r, RV):
        r_arr = np.linspace(r.ppf(0.001), r.ppf(0.999), 300)
        Er = np.trapz(r_arr ** 2 * r.pdf(r_arr), r_arr)
    else:
        Er = r ** 2
    total = total.mu_q_arr / Er
    plt.plot(w_arr, total, lw=2, color='black')

    cb = CBResidual(include_pullout=False)
    broken = SPIRRID(q=cb,
                sampling_type='PGrid',
                evars=dict(w=w_arr),
                tvars=dict(tau=tau, E_f=E_f, V_f=V_f, r=r,
                           m=m, sV0=sV0, Pf=Pf),
                n_int=1000)
    if isinstance(r, RV):
        r_arr = np.linspace(r.ppf(0.001), r.ppf(0.999), 300)
        Er = np.trapz(r_arr ** 2 * r.pdf(r_arr), r_arr)
    else:
        Er = r ** 2
    broken = broken.mu_q_arr / Er    
    plt.plot(w_arr, broken, lw=2, ls='dashed', color='black')
    plt.plot(w_arr, total-broken, lw=2, ls='dashed', color='black')   
    plt.ylim(0, 16)
    plt.xlim(0, 1.25)
    plt.xlabel('w [mm]')
    plt.ylabel('sigma_c [MPa]')
    plt.show()

from math import pi
from scipy.optimize import fsolve
def rand_xi():
    E_f, V_f, r, tau = 200e3, 0.01, 0.01, .1
    w_arr = np.linspace(0, 2.2, 300)
    T = 2. * tau / r
    for mi in [3., 7., 100.]:
        mu_xi = 0.008
        s0 = mu_xi / gamma(1. + 1./(1. + mi))
        k = np.sqrt(T/E_f)
        ef0 = k*np.sqrt(w_arr)
        G = 1 - np.exp(-(ef0/s0)**(mi+1))
        mu_int = ef0 * E_f * V_f * (1-G)
        I = s0 * gamma(1 + 1./(mi+1)) * gammainc(1 + 1./(mi+1), (ef0/s0)**(mi+1))
        mu_broken = E_f * V_f * I / (mi+1)
        plt.plot(w_arr, mu_int + mu_broken, lw = 2, color = 'black')
        #plt.plot(ef0, G, lw = 2, color = 'black')
#        wstar = (s0**(mi+1)/k**(mi+1)/mi)**(2./(mi+1))
#        numerical
#        Pf = RV('uniform', loc=0.0, scale=1.0)
#        cb = CBResidual(include_pullout=True)
#        if isinstance(r, RV):
#            r_arr = np.linspace(r.ppf(0.001), r.ppf(0.999), 300)
#            Er = np.trapz(r_arr ** 2 * r.pdf(r_arr), r_arr)
#        else:
#            Er = r ** 2
#        total = SPIRRID(q=cb,
#                sampling_type='PGrid',
#                evars=dict(w=w_arr),
#                tvars=dict(tau=tau, E_f=E_f, V_f=V_f, r=r,
#                           m=mi, sV0=sV0, Pf=Pf),
#                n_int=1000)
#        result = total.mu_q_arr / Er
#        plt.plot(w_arr, result, lw=3, color='red', ls = 'dashed', label='numerical')
    plt.xlabel('w [mm]')
    plt.ylabel('sigma_c [MPa]')
    plt.xlim(0, 2.1)
    plt.ylim(0,16)
    plt.show()

from scipy.special import gammainc, gamma

def deterministic_r():
    E_f, V_f, tau, m, sV0 = 200e3, 0.01, .1, 7.0, 3.e-3
    for ri in [0.005, 0.01, 0.015]:
        T = 2. * tau / ri
        w_arr = np.linspace(0,1.2,1000)
#         analytical solution
        s0 = ((T * (m+1) * sV0**m)/(2. * E_f * pi * ri ** 2))**(1./(m+1))
        k = np.sqrt(T/E_f)
        ef0 = k*np.sqrt(w_arr)
        G = 1 - np.exp(-(ef0/s0)**(m+1))
        mu_int = ef0 * E_f * V_f * (1-G)
        I = s0 * gamma(1 + 1./(m+1)) * gammainc(1 + 1./(m+1), (ef0/s0)**(m+1))
        mu_broken = E_f * V_f * I / (m+1)
        plt.plot(w_arr, mu_int + mu_broken, lw = 2, color = 'black')
    plt.ylim(0, 16)
    plt.xlim(0, 1.25)
    plt.show()

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax2 = ax1.twinx()

    for mi in [4.0, 5.0, 7.0]:
        r_arr = np.linspace(0.001, 0.01, 100)
        T = 2. * tau / r_arr
        s0 = ((T * (mi+1) * sV0**mi)/(2. * E_f * pi * r_arr ** 2))**(1./(mi+1))
        wstar = E_f / T * s0**2 * mi**(-2./(mi+1))
        zeta = mi**(-1./(mi+1.)) * np.exp(-1./mi) + 1. / (mi + 1.) * gamma(1 + 1./(mi+1)) * gammainc(1 + 1./(mi+1), 1./mi)
        strength = E_f * V_f * s0 * zeta
        ax1.loglog(r_arr, strength, lw=2, color='black')
        ax1.set_ylim(10,100)
        ax2.loglog(r_arr, wstar, lw=2, color='black', ls='dashed')
        ax2.set_ylim(0.1,1)
    plt.xlim(0,0.01)
    plt.show()

def rand_r():
    E_f, V_f, tau, m, sV0 = 200e3, 0.01, .1, 7.0, 3.e-3
    Pf = RV('uniform', loc=0.0, scale=1.0)
    w_arr = np.linspace(0, 1.2, 200)
    cb = CBResidual(include_pullout=True)
    for i, ri in enumerate([RV('uniform', loc=0.00999, scale=.00002), #COV = 0.01
                            RV('uniform', loc=0.005, scale=.01), #COV = 0.3
                            RV('uniform', loc=0.001, scale=.018)]): #COV = 0.5
        total = SPIRRID(q=cb,
                sampling_type='MCS',
                evars=dict(w=w_arr),
                tvars=dict(tau=tau, E_f=E_f, V_f=V_f, r=ri,
                           m=m, sV0=sV0, Pf=Pf),
                n_int=200)
        if isinstance(ri, RV):
            r_arr = np.linspace(ri.ppf(0.001), ri.ppf(0.999), 300)
            Er = np.trapz(r_arr ** 2 * ri.pdf(r_arr), r_arr)
        else:
            Er = ri ** 2
        result = total.mu_q_arr / Er
        plt.plot(w_arr, result, lw=2, color='black')
    plt.xlabel('w [mm]')
    plt.ylabel('sigma_c [MPa]')
    plt.ylim(0, 16)
    plt.xlim(0, 1.25)
    plt.show()

def deterministic_tau():
    E_f, V_f, r, m, sV0 = 200e3, 0.01, 0.01, 7., 3.e-3
    w_arr = np.linspace(0, 1.2, 300)
    for taui in [0.05, .1, .2]:
        T = 2. * taui / r
        s0 = ((T * (m+1) * sV0**m)/(2. * E_f * pi * r ** 2))**(1./(m+1))
        k = np.sqrt(T/E_f)
        ef0 = k*np.sqrt(w_arr)
        G = 1 - np.exp(-(ef0/s0)**(m+1))
        mu_int = ef0 * E_f * V_f * (1-G)
        I = s0 * gamma(1 + 1./(m+1)) * gammainc(1 + 1./(m+1), (ef0/s0)**(m+1))
        mu_broken = E_f * V_f * I / (m+1)
        plt.plot(w_arr, mu_int + mu_broken, lw = 2, color = 'black')
    plt.ylim(0, 16)
    plt.xlim(0, 1.25)
    plt.show()

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax2 = ax1.twinx()

    for mi in [4.0, 5.0, 7.0]:
        tau_arr = np.linspace(0.05, .2, 5)
        T = 2. * tau_arr / r
        s0 = ((T * (mi+1) * sV0**mi)/(2. * E_f * pi * r ** 2))**(1./(mi+1))
        wstar = E_f / T * s0**2 * mi**(-2./(mi+1))
        zeta = mi**(-1./(mi+1.)) * np.exp(-1./mi) + 1. / (mi + 1.) * gamma(1 + 1./(mi+1)) * gammainc(1 + 1./(mi+1), 1./mi)
        strength = E_f * V_f * s0 * zeta
        ax1.loglog(tau_arr, strength, lw=2, color='black')
        ax1.set_ylim(0,15)
        ax2.loglog(tau_arr, wstar, lw=2, color='black', ls='dashed')
        ax2.set_ylim(0.1,1.)
    plt.xlim(0.05,.2)
    plt.ylim(0)
    plt.show()

def rand_tau():
    E_f, V_f, r, m, sV0 = 200e3, 0.01, 0.01, 7.0, 3.e-3
    Pf = RV('uniform', loc=0.0, scale=1.0)
    w_arr = np.linspace(0, 1., 200)
    cb = CBResidual(include_pullout=True)
    for i, taui in enumerate([RV('uniform', loc=.0999, scale=.0002), #COV = 0.01
                            RV('uniform', loc=.05, scale=.1), #COV = 0.3
                            RV('uniform', loc=.01, scale=0.18)]): #COV = 0.5
        total = SPIRRID(q=cb,
                sampling_type='MCS',
                evars=dict(w=w_arr),
                tvars=dict(tau=taui, E_f=E_f, V_f=V_f, r=r,
                           m=m, sV0=sV0, Pf=Pf),
                n_int=300)
        if isinstance(r, RV):
            r_arr = np.linspace(r.ppf(0.001), r.ppf(0.999), 300)
            Er = np.trapz(r_arr ** 2 * r.pdf(r_arr), r_arr)
        else:
            Er = r ** 2
        result = total.mu_q_arr / Er
        plt.plot(w_arr, result, lw=2, color='black')
    plt.xlabel('w [mm]')
    plt.ylabel('sigma_c [MPa]')
    plt.ylim(0, 16)
    plt.xlim(0, 1.25)
    plt.show()

def g_ell():
    E_f, r, tau = 200e3, 0.01, .1
    w = 0.3
    T = 2. * tau / r
    e = np.sqrt(T*w/E_f)
    a = e / 2. / tau * r * E_f
    z_arr = np.linspace(0.0, a, 500)
    m0 = .7
    pdf = m0/a * (1-z_arr/a)**(m0-1)
    plt.plot(z_arr, pdf, lw=2, ls='dashed', color='black')
    m1 = 5.0
    pdf = m1/a * (1-z_arr/a)**(m1-1)
    plt.plot(z_arr, pdf, lw=2, color='black')
    m2 = 6.0
    pdf = m2/a * (1-z_arr/a)**(m2-1)
    plt.plot(z_arr, pdf, lw=2, color='black')
    m3 = 7.0
    pdf = m3/a * (1-z_arr/a)**(m3-1)
    plt.plot(z_arr, pdf, lw=2, color='black')
    plt.ylim(0, 0.15)
    #plt.xlim(0, 20.)
    plt.show()

def mu_ell():
    E_f, r, tau, sV0 = 200e3, 0.01, .1, 3.e-3
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax2 = ax1.twinx()
    def CDFa(e, mm):
        s = sV0*(pi*r**2)**(-1./mm)
        T = 2. * tau / r
        a = e / T
        return 1 - np.exp(-a * 2 * E_f * (e / s) ** mm / (mm + 1))
    w_arr = np.linspace(0.0, 1.3, 1000)
    T = 2. * tau / r
    e_arr = np.sqrt(T*w_arr/E_f)
    def mu_ell(mm):
        n = mm+1
        s0 = ((T * (mm+1) * sV0**mm)/(2. * E_f * pi * r ** 2))**(1./(mm+1))
        I = s0 * gamma(1. + 1./n) * gammainc(1. + 1./n, (e_arr/s0)**n)
        mu_ell = I * E_f / T / n / CDFa(e_arr, mm)
        ax1.plot(w_arr, mu_ell, lw=2, color='black')
    mu_ell(5.0)
    ax2.plot(w_arr, CDFa(e_arr, 5.0), lw=2, ls='dashed', color='black')
    mu_ell(6.0)
    ax2.plot(w_arr, CDFa(e_arr, 6.0), lw=2, ls='dashed', color='black')
    mu_ell(7.0)
    ax2.plot(w_arr, CDFa(e_arr, 7.0), lw=2, ls='dashed', color='black')
    plt.xlim(0,1.2)
    plt.show()

#general_diagram()
#rand_xi()
#deterministic_r()
#rand_r()
deterministic_tau()
#rand_tau()
#g_ell()
#mu_ell()
