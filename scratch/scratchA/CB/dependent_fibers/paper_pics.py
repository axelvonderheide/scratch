'''
Created on Nov 19, 2012

@author: rostar
'''

from etsproxy.mayavi import mlab as m
import numpy as np
from scipy.stats import weibull_min
from stats.spirrid import make_ogrid as orthogonalize
from composite_CB_model import CompositeCrackBridge
from composite_CB_postprocessor import CompositeCrackBridgeView
from matplotlib import pyplot as plt
from stats.spirrid.rv import RV
from reinforcement import Reinforcement, WeibullFibers
from stats.spirrid.spirrid import SPIRRID
from paper_pics.rf_rigid_mtrx import CBRigidMatrix, CBRigidMatrixSP

reinf1 = Reinforcement(r=0.00345,#RV('uniform', loc=0.001, scale=0.005),
                      tau=RV('uniform', loc=.5, scale=.2),
                      V_f=0.3,
                      E_f=70e3,
                      xi=WeibullFibers(shape=5., scale=0.02, L0=10.),#RV('weibull_min', shape=5., scale=.02),
                      n_int=15,
                      label='AR glass')

reinf2 = Reinforcement(r=0.003,#RV('uniform', loc=0.002, scale=0.002),
                      tau=0.5,#RV('uniform', loc=.5, scale=.1),
                      V_f=0.2,
                      E_f=200e3,
                      xi=RV('weibull_min', shape=5., scale=.02),
                      n_int=15,
                      label='carbon')

model = CompositeCrackBridge(E_m=25e3,
                             reinforcement_lst=[reinf2],
                             Ll=50.,
                             Lr=50.)

ccb_view = CompositeCrackBridgeView(model=model)


def ld_rigid_vs_el_mtrx():
    w_arr = np.linspace(0.0, 0.8, 200)
    r = 0.003
    tau = 0.3
    V_f = 0.11
    E_f = 200e3
    n_int = 30
    E_m = 25e3
    xi = 0.02
    Xxi = RV('weibull_min', shape=5., scale=.02)
    Xr = RV('uniform', loc=.002, scale=.002)
    Xtau = RV('uniform', loc=.02, scale=.98)
    reinf1 = Reinforcement(r=r, tau=tau, V_f=V_f, E_f=E_f, xi=xi, n_int=n_int)
    reinf2 = Reinforcement(r=r, tau=tau, V_f=V_f, E_f=E_f, xi=Xxi, n_int=100)
    reinf3 = Reinforcement(r=Xr, tau=Xtau, V_f=V_f, E_f=E_f, xi=Xxi, n_int=15)
    model1 = CompositeCrackBridge(E_m=E_m, reinforcement_lst=[reinf1], Ll=35., Lr=35.)
    model2 = CompositeCrackBridge(E_m=E_m, reinforcement_lst=[reinf2], Ll=35., Lr=35.)
    model3 = CompositeCrackBridge(E_m=E_m, reinforcement_lst=[reinf3], Ll=35., Lr=35.)
    ccb_view1 = CompositeCrackBridgeView(model=model1)
    ccb_view2 = CompositeCrackBridgeView(model=model2)
    ccb_view3 = CompositeCrackBridgeView(model=model3)
    sigma_c1, u1, sigma_c2, u2, sigma_c3, u3 = [], [], [], [], [], []
    for w in w_arr:
        ccb_view1.model.w = w/2
        ccb_view2.model.w = w
        ccb_view3.model.w = w
        u1.append(ccb_view1.u_evaluated)
        sigma_c1.append(ccb_view1.sigma_c)
        u2.append(ccb_view2.u_evaluated)
        sigma_c2.append(ccb_view2.sigma_c)
        u3.append(ccb_view3.u_evaluated)
        sigma_c3.append(ccb_view3.sigma_c)
    plt.plot(w_arr / 2., sigma_c1, lw = 2, color='black', label='w')
    plt.plot(w_arr, sigma_c2, lw = 2, color='black')
    plt.plot(w_arr, sigma_c3, lw = 2, color='black')
    plt.plot(u1, sigma_c1, ls='dashed', lw = 2, color='black', label='u')
    plt.plot(u2, sigma_c2, ls='dashed', lw = 2, color='black')
    plt.plot(u3, sigma_c3, ls='dashed', lw = 2, color='black')
    plt.xlabel('w, u [mm]')
    plt.ylabel('$\sigma_c$ [MPa]')

def profiles_rigid_vs_el_mtrx():
    r = 0.003
    tau = 0.3
    V_f = 0.11
    E_f = 200e3
    E_m = 25e10
    Ll = 35.
    Lr = 35.
    Xxi = RV('weibull_min', shape=5., scale=.02)
    Xr = RV('uniform', loc=.002, scale=.002)
    Xtau = RV('uniform', loc=.02, scale=.98)
    reinf2 = Reinforcement(r=Xr, tau=Xtau, V_f=V_f, E_f=E_f, xi=Xxi, n_int=15)
    model2 = CompositeCrackBridge(E_m=E_m, reinforcement_lst=[reinf2], Ll=Ll, Lr=Lr)
    ccb_view2 = CompositeCrackBridgeView(model=model2)

    ccb_view2.sigma_c_max
    x2 = np.hstack((-Ll, ccb_view2.x_arr, Lr))
    mu_epsf2 = np.hstack((ccb_view2.mu_epsf_arr[0], ccb_view2.mu_epsf_arr, ccb_view2.mu_epsf_arr[-1]))
    epsm2 = np.hstack((ccb_view2.epsm_arr[0], ccb_view2.epsm_arr, ccb_view2.epsm_arr[-1]))
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(x2, mu_epsf2, lw=2, color='black')
    ax.plot(x2, epsm2, lw=2, color='black')
    p = ax.fill_between(x2, epsm2, mu_epsf2, facecolor='none')
    from matplotlib.patches import PathPatch
    path = p.get_paths()[0]
    p1 = PathPatch(path, fc="none", hatch="/")
    ax.add_patch(p1)
    p1.set_zorder(p.get_zorder()-0.1)
    p = ax.fill_between(x2, mu_epsf2, facecolor='none')
    path = p.get_paths()[0]
    p1 = PathPatch(path, fc="none", hatch="\\")
    ax.add_patch(p1)
    p1.set_zorder(p.get_zorder()-0.1)
    
    reinf1 = Reinforcement(r=r, tau=tau, V_f=V_f, E_f=E_f, xi=Xxi, n_int=15)
    model1 = CompositeCrackBridge(E_m=E_m, reinforcement_lst=[reinf1], Ll=Ll, Lr=Lr)
    ccb_view1 = CompositeCrackBridgeView(model=model1)
    ccb_view1.sigma_c_max
    x1 = np.hstack((-Ll, ccb_view1.x_arr, Lr))
    mu_epsf1 = np.hstack((ccb_view1.mu_epsf_arr[0], ccb_view1.mu_epsf_arr, ccb_view1.mu_epsf_arr[-1]))
    epsm1 = np.hstack((ccb_view1.epsm_arr[0], ccb_view1.epsm_arr, ccb_view1.epsm_arr[-1]))
    ax.plot(x1, mu_epsf1, lw=2, color='black')
    ax.plot(x1, epsm1, lw=2, color='black')
    
    plt.xlabel('z [mm]')
    plt.xlim(-35, 35)
    plt.ylim(0)
    plt.ylabel('$\epsilon$ [-]')

def random_domain(w):
    Ef = 70e3
    Fxi = weibull_min(5., scale = 0.02)
    r = np.linspace(0.001, 0.005, 100)
    tau = np.linspace(0., 1., 100)
    e_arr = orthogonalize([np.arange(len(r)), np.arange(len(tau))])

    tau = tau.reshape(1, len(tau))
    r = r.reshape(len(r), 1)
    eps0 = np.sqrt(w * 2 * tau / r / Ef)
    F = Fxi.cdf(eps0)
    a = r * Ef / 2. / tau
    m.surf(e_arr[0], e_arr[1], 50 * F / np.max(F))
    m.surf(e_arr[0], e_arr[1], 50 * eps0 / np.max(eps0))
    m.surf(e_arr[0], e_arr[1], np.nan_to_num(a) / 100.)
    m.show()

def sigma_f(w_arr):
    sf_arr = ccb_view.sigma_f_lst(w_arr)
    #for i, reinf in enumerate(ccb_view.model.reinforcement_lst):
    #    plt.plot(w_arr, sf_arr[:, i], label=reinf.label)
    plt.plot(w_arr, sf_arr[:,0], ls='dashed', label='$\sigma_{\mathrm{f}1}$', color='black')
    plt.plot(w_arr, sf_arr[:,1], ls='dotted', label='$\sigma_{\mathrm{f}2}$', color='black')
    plt.xlabel('crack opening w [mm]')
    plt.ylabel('stress [MPa]')
    plt.ylim(0)

def sigma_c_w(w_arr, r, tau, E_f, E_m, V_f, xi, n_int):
    reinf = Reinforcement(r=r, tau=tau, E_f=E_f, V_f=V_f,
                          xi=xi, n_int=n_int)
    model = CompositeCrackBridge(E_m=E_m,
                             reinforcement_lst=[reinf],
                             Ll=1000.,
                             Lr=1000.)
    ccb_view = CompositeCrackBridgeView(model=model)
    sigma_c = []
    for w in w_arr:
        ccb_view.model.w = w
        sigma_c.append(ccb_view.sigma_c)
    plt.plot(w_arr, sigma_c, lw=2, label='elastic_mtrx')

    cb = CBRigidMatrix()
    spirrid = SPIRRID(q=cb, sampling_type='PGrid',
                      evars=dict(w=w_arr),
                      tvars=dict(r=r, tau=tau, E_f=E_f, V_f=V_f,
                                 xi=xi),
                      n_int=n_int)
    plt.plot(w_arr, E_f * V_f * spirrid.mu_q_arr, lw=2, label='rigid_mtrx')

def errors(k_ratio):
    sigmax_el = []
    wmax = []
    for k in k_ratio:
        E_m = E_f * V_f / k / (1 - V_f)
        xi = WeibullFibers(shape=5., scale=0.02, L0=10.)#RV('weibull_min', shape=50., scale=.02)
        r = 0.002#RV('uniform', loc=0.002, scale=.002)
        tau = RV('weibull_min', shape=80., scale=.5)
        reinf = Reinforcement(r=r, tau=tau, V_f=0.1, E_f=200e3, xi=xi, n_int=20, label='carbon')
        model = CompositeCrackBridge(E_m=E_m,
                                 reinforcement_lst=[reinf],
                                 Ll=100.,
                                 Lr=100.)
        ccb_view = CompositeCrackBridgeView(model=model)
        sigmax_el.append(ccb_view.sigma_c_max[0])
        wmax.append(ccb_view.sigma_c_max[1])
    ccb_view.model.E_m = 10e10
    rigid_sigma = ccb_view.sigma_c_max[0]
    rigid_w = ccb_view.sigma_c_max[1]

    plt.plot(k_ratio, np.array(sigmax_el)/rigid_sigma, label='$\sigma_\mathrm{c, max, el}/\sigma_\mathrm{c, max, rigid}$', lw=2, color='red')
    plt.ylabel('$\sigma_\mathrm{c, max, el}/\sigma_\mathrm{c, max, rigid}$ and $w_\mathrm{max, el} / w_\mathrm{max, rigid}$')
    plt.xlabel('$K_\mathrm{f}/K_\mathrm{m}$')
    plt.title('effect of matrix stiffness')
    plt.plot(k_ratio, np.array(wmax)/rigid_w, label='$w_\mathrm{max, el} / w_\mathrm{max, rigid}$', lw=2, color='blue')
    plt.ylim(0, 1.3)

def profile(wr, we, r, tau, E_f, E_m, V_f, xi, n_int):
    reinf = Reinforcement(r=r, tau=tau, E_f=E_f, V_f=V_f,
                          xi=xi, n_int=n_int)
    model = CompositeCrackBridge(w=we,
                                 E_m=E_m,
                                 reinforcement_lst=[reinf],
                                 Ll=1000., Lr=1000.)
    ccb_view = CompositeCrackBridgeView(model=model)
    plt.plot(ccb_view.x_arr, ccb_view.mu_epsf_arr, lw=2, label='el_mtrx_yarn')
    plt.plot(ccb_view.x_arr, ccb_view.epsm_arr, lw=2, label='el_mtrx_mtrx')

    cb = CBRigidMatrixSP()
    x_arr = np.hstack((-np.linspace(3.*ccb_view.x_arr[-1], 0., 100),
                        np.linspace(0.0, 3.*ccb_view.x_arr[-1],100)))
    spirrid = SPIRRID(q=cb, sampling_type='PGrid',
                      evars=dict(x=x_arr),
                      tvars=dict(w=wr, r=r, tau=tau, E_f=E_f, V_f=V_f, xi=xi),
                      n_int=n_int)
    plt.plot(x_arr, spirrid.mu_q_arr, lw=2, label='rigid_mtrx_yarn')
    plt.xlabel('position [mm]')
    plt.ylabel('strain')

def sigma_c_u(w_arr, r, tau, E_f, E_m, V_f, xi, n_int):
    reinf = Reinforcement(r=r, tau=tau, E_f=E_f, V_f=V_f,
                          xi=xi, n_int=n_int)
    model = CompositeCrackBridge(E_m=E_m,
                             reinforcement_lst=[reinf],
                             Ll=30.,
                             Lr=30.)
    ccb_view = CompositeCrackBridgeView(model=model)
    sigma_c = []
    u = []
    for w in w_arr:
        ccb_view.model.w = w
        u.append(ccb_view.u_evaluated)
        sigma_c.append(ccb_view.sigma_c)
    plt.plot(w_arr, sigma_c, lw=1, color='blue')
    plt.plot(u, sigma_c, lw=1, color='red')

    cb = CBRigidMatrix()
    spirrid = SPIRRID(q=cb, sampling_type='PGrid',
                      evars=dict(w=w_arr),
                      tvars=dict(r=r, tau=tau, E_f=E_f, V_f=V_f,
                                 xi=xi),
                      n_int=n_int)
    plt.plot(w_arr, E_f * V_f * spirrid.mu_q_arr, lw=2, color='black')

def integ_domain(isolines):
    #T = np.linspace(2., 250., isolines)
    tau = np.linspace(0.01, 0.5)
    #for t in T:
    plt.plot(tau, 2 * tau / 100., lw=2., color='black')
    plt.xlabel('')
    plt.ylabel('r')

#ld_rigid_vs_el_mtrx()
#profiles_rigid_vs_el_mtrx()
#elastic_matrix(np.linspace(.0, .5, 100))
#sigma_f(np.linspace(.0, .3, 100))
#rigid_mtrx()
#plt.subplot(1,3,1)
#errors(np.linspace(0.000001, 2., 20))
#plt.legend(loc='best')
#plt.subplot(1,3,2)
#sigma_c_w(w_arr=np.linspace(.0, .3, 100),
#          r=0.002, tau=RV('weibull_min', shape=5., scale=.5),
#          E_f=200e3, E_m=25e3,
#          V_f=0.9, xi=RV('weibull_min', shape=5., scale=.02),
#          n_int=50)
#plt.legend(loc='best')
#plt.subplot(1,3,3)
#profile(wr=.087, we=0.067, r=0.002, tau=RV('weibull_min', shape=5., scale=.5),
#        E_f=200e3, E_m=25e3,
#        V_f=0.05, xi=RV('weibull_min', shape=5., scale=.02),
#        n_int=50)
#for em in np.linspace(5e3, 100e3, 4):
#    sigma_c_u(w_arr=np.linspace(.0, .5, 100),
#          r=0.002, tau=RV('weibull_min', shape=5., scale=.5),
#          E_f=200e3, E_m=em,
#          V_f=0.1, xi=RV('weibull_min', shape=5., scale=.02),
#          n_int=50)
#integ_domain(10)
plt.legend(loc='best')
plt.show()
#random_domain(0.15)