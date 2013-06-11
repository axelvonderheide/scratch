#---#-------------------------------------------------------------------------------
#
# Copyright (c) 2009, IMB, RWTH Aachen.
# All rights reserved.
#
# This software is provided without warranty under the terms of the BSD
# license included in simvisage/LICENSE.txt and may be redistributed only
# under the conditions described in the aforementioned license.  The license
# is also available online at http://www.simvisage.com/licenses/BSD.txt
#
# Thanks for using Simvisage open source!
#
# Created on Jun 14, 2010 by: rch

import numpy as np

from matplotlib import pyplot as plt
from scipy.optimize import brentq


def H(x):
    return x >= 0.0


class WAnalytIter():

    def crackbridge(self, w, l, T, Kf, Km, Vf, omega):
        #Phase A : Both sides debonding .
        Kc = Kf + Km
        c = Kc * T * l / 2.
        q0 = (np.sqrt(c ** 2 + w * Kf * Km * Kc * T) - c) / Km
        return q0 / Vf / (1. - omega)

    def pullout(self, w, l, T, Kf, Km, Vf, Lmin, Lmax, omega):
        #Phase B : Debonding of shorter side is finished
        Kc = Kf + Km
        c = Kc * T * (Lmin + l)
        f = T ** 2 * Lmin ** 2 * Kc ** 2
        q1 = (np.sqrt(c ** 2. + f + 2 * w * Kc * T * Kf * Km) - c) / Km
        return q1 / Vf / (1. - omega)

    def linel(self, w, l, T, Kf, Km, Vf, Lmax, Lmin, omega):
        #Phase C: Both sides debonded - linear elastic behavior.
        Kc = Kf + Km
        q2 = (2. * w * Kf * Km + T * Kc * (Lmin ** 2 + Lmax ** 2))/(2. * Km * (Lmax + l + Lmin))
        return q2 / Vf / (1. - omega)

    def q(self, w, tau, l, E_f, E_m, theta, xi, phi, Ll, Lr, V_f, r, omega):
        #assigning short and long embedded length
        Lmin = np.minimum(Ll, Lr)
        Lmax = np.maximum(Ll, Lr)

        Lmin = np.maximum(Lmin - l / 2., 1e-15)
        Lmax = np.maximum(Lmax - l / 2., 1e-15)

        #maximum condition for free length
        l = np.minimum(l / 2., Lr) + np.minimum(l / 2., Ll)

        #defining variables
        w = w - theta * l
        l = l * (1 + theta)
        w = H(w) * w
        T = 2. * tau * V_f * (1. - omega) / r
        Km = (1. - V_f) * E_m
        Kf = V_f * (1. - omega) * E_f
        Kc = Km + Kf

        # double sided debonding
        #q0 = self.crackbridge(w, l, T, Kr, Km, Lmin, Lmax)
        q0 = self.crackbridge(w, l, T, Kf, Km, V_f, omega)

        # displacement at which the debonding to the closer clamp is finished
        w0 = (Lmin+l) * Lmin * Kc * T / Kf / Km

        # debonding of one side; the other side is clamped
        q1 = self.pullout(w, l, T, Kf, Km, V_f, Lmin, Lmax, omega)

        # displacement at which the debonding is finished at both sides
        e1 = Lmax * Kc * T / Km / Kf
        w1 = e1 * (l + Lmax / 2.) + (e1 + e1 * Lmin / Lmax) * Lmin / 2.

        # debonding completed at both sides, response becomes linear elastic
        q2 = self.linel(w, l, T, Kf, Km, V_f, Lmax, Lmin, omega)

        # cut out definition ranges
        q0 = H(w) * (q0 + 1e-15) * H(w0 - w)
        q1 = H(w - w0) * q1 * H(w1 - w)
        q2 = H(w - w1) * q2

        #add all parts
        q = q0 + q1 + q2
        return q
    
    def damage_func(self, w, tau, l, E_f, E_m, theta, xi, phi, Ll, Lr, V_f, r, omega):
        # include breaking strain
        q = self.q(w, tau, l, E_f, E_m, theta, xi, phi, Ll, Lr, V_f, r, omega)
        damage = lambda x: xi._distr.cdf(x)
        return damage(q / E_f)

    def residuum(self, w, tau, l, E_f, E_m, theta, xi, phi, Ll, Lr, V_f, r, omega):
        damage_func = self.damage_func(w, tau, l, E_f, E_m, theta, xi, phi, Ll, Lr, V_f, r, omega)
        return damage_func - omega - 1e-15

    def eval_w_q(self, tau, l, E_f, E_m, theta, xi, phi, Ll, Lr, V_f, r, ctrl_damage):
        w_lst = []
        q_lst = []
        for omega in ctrl_damage:
            wD = brentq(self.residuum, 0.0, 5.0,
                        args=(tau, l, E_f, E_m, theta, xi, phi, Ll, Lr, V_f, r, omega))
            q = self.q(wD, tau, l, E_f, E_m,
                                 theta, xi, phi, Ll, Lr, V_f, r, omega) * (1. - omega)
            w_lst.append(wD)
            q_lst.append(q)
        return np.array(w_lst), np.array(q_lst)

if __name__ == '__main__':
    from stats.spirrid.rv import RV
    a = WAnalytIter()
    # filaments
    r = 0.00345
    V_f = 0.04
    tau = 0.5#RV('uniform', loc=0.2, scale=1.)
    E_f = 200e3
    E_m = 25e3
    l = 0.0
    theta = 0.0
    phi = 1.
    Ll = 100.
    Lr = 100.
    xi = RV('weibull_min', shape=5., scale=.02)
    ctrl_damage = np.linspace(0.0, .9999, 50)
    w = np.linspace(0.,0.7,100)
    #q = a.damage_func(w, tau, l, E_f, E_m, theta, xi, phi, Ll, Lr, V_f, r, 0.3)
    w, q = a.eval_w_q(tau, l, E_f, E_m, theta, xi, phi, Ll, Lr, V_f, r, ctrl_damage)
    plt.plot(w, q)
    plt.show()
