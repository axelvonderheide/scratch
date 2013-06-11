'''
Created on Oct 26, 2012

@author: rostar
'''

from composite_CB_model import CompositeCrackBridge
from composite_CB_modelview import CompositeCrackBridgeView
from reinforcement import Reinforcement
from stats.spirrid.rv import RV
from matplotlib import pyplot as plt
import numpy as np

if __name__ == '__main__':

    # AR-glass
    reinf1 = Reinforcement(r=RV('uniform', loc=0.012, scale=0.002),
                          tau=RV('uniform', loc=.3, scale=.1),
                          V_f=0.2,
                          E_f=72e3,
                          xi=RV('weibull_min', shape=20., scale=.02),
                          n_int=15)

    # carbon
    reinf2 = Reinforcement(r=RV('uniform', loc=0.002, scale=0.002),
                          tau=RV('uniform', loc=.6, scale=.1),
                          V_f=0.05,
                          E_f=200e3,
                          xi=RV('weibull_min', shape=10., scale=.015),
                          n_int=15)

    # instance of CompCrackBridge with matrix E and BC
    model = CompositeCrackBridge(E_m=25e3,
                               reinforcement_lst=[reinf1, reinf2],
                               Ll=10.,
                               Lr=20.)

    ccb_view = CompositeCrackBridgeView(model=model)

    def profile(w):
        '''evaluates and plots the strain profile in matrix
        and the mean strain profile in the reinforcement'''
        ccb_view.model.w = w
        plt.plot(ccb_view.x_arr, ccb_view.em_arr, label='matrix')
        plt.plot(ccb_view.x_arr, ccb_view.mu_epsf_arr, label='yarn')
        plt.xlabel('position [mm]')
        plt.ylabel('strain [-]')
        plt.legend(loc='best')

    def norm_stress_w(w_arr):
        '''evaluates the normalized stress
        (force acting on unit CS) vs crack opening'''
        norm_stress = []
        w_err = []
        for w in w_arr:
            ccb_view.model.w = w
            print 'progress: ', 100 * w / w_arr[-1], '%'
            norm_stress.append(ccb_view.max_norm_stress)
            w_err.append((ccb.w_evaluated - ccb.w) / (ccb.w + 1e-10))
        plt.figure()
        plt.plot(w_arr, w_err, label='error in w')
        plt.xlabel('control w')
        plt.ylabel('rel. error in w [$w_{eval}/w_{ctrl}$]')
        plt.legend(loc='best')
        plt.figure()
        plt.plot(w_arr, norm_stress, lw=2, label='l-d')
        plt.xlabel('w [mm]')
        plt.ylabel('norm stress [MPa]')
        plt.legend(loc='best')

    def bundle_comparison(w_arr, L, shape, scale, E):
        '''bundle (Weibull fibers) response for comparison with the CB model'''
        from scipy.stats import weibull_min
        eps = w_arr / L * (1. - weibull_min(shape, scale=scale).cdf(w_arr / L))
        plt.plot(w_arr / L, eps * E, lw=4, color='red', ls='dashed', label='FB model')

        bundle = Reinforcement(r=0.13, tau=0.00001, V_f=0.9999, E_f=E,
                          xi=RV('weibull_min', shape=shape, scale=scale),
                          n_int=50)
        ccb = CompositeCrackBridge(E_m=25e3,
                                   reinforcement_lst=[bundle],
                                   Ll=L / 2.,
                                   Lr=L / 2.)
        eps = []
        for w in w_arr:
            ccb.w = w
            eps.append(ccb.max_norm_stress / E)
        plt.plot(w_arr / L, np.array(eps) * E, color='blue', lw=2, label='CB model')
        plt.legend(loc='best')

    def analytical_comparison():
        '''for the case when tau is deterministic,
        there is an analytical solution.
        The differences are caused by the additional matrix
        stiffness due to broken fibers, which are in the CB model
        added to matrix stiffness. As the matrix E grows and the V_f
        decreases, the solutions tend to get closer'''
        tau, E_f, E_m, V_f = 0.2, 200e3, 25e3, 0.3
        r, shape, scale = 0.00345, 5., 0.02

        # analytical solution for damage controlled test
        ctrl_damage = np.linspace(0.0, .99, 100)

        def crackbridge(w, tau, E_f, E_m, V_f, r, omega):
            Kf = E_f * V_f * (1 - omega)
            Km = E_m * (1 - V_f) + E_f * V_f * omega 
            Kc = Kf + Km
            T = 2. * tau * V_f * (1. - omega) / r
            c = np.sqrt(Kc * T / Km / Kf)
            return c * np.sqrt(w) * (1 - omega)

        def w_omega(tau, E_f, E_m, V_f, r, omega, shape, scale):
            Kf = E_f * V_f * (1 - omega)
            Km = E_m * (1 - V_f) + E_f * V_f * omega
            Kc = Kf + Km
            T = 2. * tau * V_f * (1. - omega) / r
            return (-np.log(1. - omega)) ** (2. / shape) \
                    * scale ** 2 * Km * Kf / Kc / T

        w_lst = [w_omega(tau, E_f, E_m, V_f, r, omega, shape, scale)
                 for omega in ctrl_damage]
        epsf = crackbridge(np.array(w_lst), tau, E_f, E_m, V_f, r, ctrl_damage)
        plt.plot(np.array(w_lst), epsf * E_f, color='red',
                 lw=4, ls='dashed', label='analytical')

        reinf = Reinforcement(r=r, tau=tau, V_f=V_f, E_f=E_f,
                              xi=RV('weibull_min', shape=shape, scale=scale),
                              n_int=20)

        ccb = CompositeCrackBridge(E_m=E_m, reinforcement_lst=[reinf],
                                   Ll=1000., Lr=1000.)

        stress = []
        w_arr = np.linspace(0.0, np.max(w_lst), 100)
        for w in w_arr:
            ccb.w = w
            stress.append(ccb.max_norm_stress)
        plt.plot(w_arr, stress, color='blue', lw=2, label='CB model')
        plt.legend(loc='best')

    #profile(.01)
    #norm_stress_w(np.linspace(.0, .38, 50))
    #bundle_comparison(np.linspace(0, 0.65, 30), 20., 5., 0.02, 70e3)
    analytical_comparison()
    plt.show()
