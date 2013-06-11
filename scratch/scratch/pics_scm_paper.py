'''
Created on 24.1.2012

@author: Q
'''
from quaducom.micro.resp_func.cb_emtrx_clamped_fiber import CBEMClampedFiberSP, CBEMClampedFiber, H
import numpy as np
from matplotlib import pyplot as plt
from math import pi
from scipy.stats import uniform, norm, weibull_min
from stats.spirrid.spirrid import SPIRRID
from stats.spirrid.rv import RV
from stats.spirrid import make_ogrid
from etsproxy.mayavi import mlab as m

tau = .2
Af = 5.31e-4
Ef = 72e3
Am = 50. / 1700.
Em = 30e3
l = 0.
theta = 0.0
xi = 5000.5
phi = 1.
Ll = 50.
Lr = 50.
Nf = 1.

def Pw():
    w = np.linspace(0, .5, 300)

    def crackbridge(w, l, T, Kr, Km, Lmin, Lmax):
        c = (2 * l * (Km + Kr) + Kr * (Lmin + Lmax))
        P0 = (T * (Km + Kr)) / (2. * Km ** 2) * (np.sqrt(c ** 2 + 4 * w * Kr * Km ** 2 / T) - c)
        return P0

    def pullout(u, l, T, Kr, Km, L):
        c = l * (Km + Kr) + L * Kr
        P1 = (T * (Km + Kr)) / Km ** 2 * (np.sqrt(c ** 2 + 2 * u * Kr * Km ** 2 / T) - c)
        return P1

    def linel(u, l, T, Kr, Km, L):
        P2 = (T * L ** 2 + 2 * u * Kr) / 2. / (L + l)
        return P2

    def recall(w, tau, l, A_f, E_f, E_m, A_m, theta, xi, phi, Ll, Lr, Nf):
        Lmin = np.minimum(Ll, Lr)
        Lmax = np.maximum(Ll, Lr)
        Lmin = np.maximum(Lmin - l / 2., 0)
        Lmax = np.maximum(Lmax - l / 2., 0)
        l = np.minimum(l / 2., Lr) + np.minimum(l / 2., Ll)
        l = l * (1 + theta)
        w = w - theta * l
        w = H(w) * w
        o = np.sqrt(4. * A_f * Nf * pi)
        T = tau * phi * o
        Km = A_m * E_m
        Kr = A_f * E_f
        l0 = l / 2.
        q0 = crackbridge(w, l0, T, Kr, Km, Lmin, Lmax)
        w0 = T * Lmin * ((2 * l0 + Lmin) * (Kr + Km) + Kr * Lmax) / (Km * Kr)
        Q0 = Lmin * T * (Km + Kr) / (Km)
        l1 = 2 * Lmin + l
        L1 = Lmax - Lmin
        q1 = pullout(w - w0, l1, T, Kr, Km, L1) + Q0
        w1 = (L1 * T * (Kr + Km) * (L1 + l1)) / Kr / Km - T * L1 ** 2 / 2. / Kr
        q2 = linel(w - w0, l1, T, Kr, Km, L1) + Q0
        q0 = H(w) * (q0 + 1e-15) * H(w0 - w)
        q1 = H(w - w0) * q1 * H(w1 + w0 - w)
        q2 = H(w - w1 - w0) * q2
        q = q0 + q1 + q2
        q = q * H(Kr * xi - q)
        mask0 = q0 * H(Kr * xi - q0) > 0
        mask1 = q1 * H(Kr * xi - q1) > 0
        mask2 = q2 > 0
        return (mask0, q0[mask0]), (mask1, q1[mask1]), (mask2, (q2 * H(Kr * xi - q2))[mask2])
    q = recall(w, tau, l, Af, Ef, Em, Am, theta, xi, phi, Ll, Lr, Nf)
    plt.plot(w[q[0][0]], q[0][1], lw = 2, ls = '-',
              color = 'black', label = 'stage A')
    plt.plot(w[q[1][0]], q[1][1], lw = 2, ls = '--',
              color = 'black', label = 'stage B')
    plt.plot(w[q[2][0]], q[2][1], lw = 2, ls = 'dotted',
              color = 'black', label = 'stage C')
    plt.xlabel('$\mathrm{crack width} \, w \mathrm{[mm]}$', fontsize = 16)
    plt.ylabel('$\mathrm{force} \, P_\mathrm{f,0} \mathrm{[N]}$', fontsize = 16)
    plt.title('$\mathrm{filament \, crack \, bridge}$' , fontsize = 20)
    plt.legend(loc = 'best')
    plt.show()

def SP():
    cbcsp = CBEMClampedFiberSP()
    x = np.linspace(-50, 20, 100)
    q = cbcsp(.5, x, tau, l, Af, Ef, Em, Am, theta, xi, phi, Ll, Lr, Nf)
    plt.plot(x, q, lw = 4, color = 'black', label = 'force along filament')
    plt.xticks(fontsize = 20)
    plt.yticks(fontsize = 20)
    #plt.legend( loc = 'best' )
    plt.xlabel('$\mathrm{position} \, x \mathrm{[mm]}$', fontsize = 30)
    plt.ylabel('$\mathrm{force} \, P_\mathrm{f} \mathrm{[N]}$', fontsize = 30)
    plt.ylim(0, 0.5)
    plt.show()

def random_samples_Pw(n):
    w = np.linspace(0, 1.7, 300)
    P = CBEMClampedFiber()
    spirrid = SPIRRID(q = P,
                      sampling_type = 'LHS',
                      evars = dict(w = w
                                    ),
                      tvars = dict(Ll = 50.,
                                   Lr = 20.,
                                   tau = RV('uniform', 0.05, .15),
                                   l = RV('uniform', 5.0, 10.0),
                                   A_f = Af,
                                   E_f = Ef,
                                   theta = RV('uniform', 0.0, .02),
                                   xi = RV('weibull_min', scale = 0.017, shape = 5, n_int = 10),
                                   phi = phi,
                                   E_m = Em,
                                   A_m = 50.,
                                   Nf = 1700.
                                    ),
                    n_int = 20)

    for i in range(n):
        if i == n - 1:
            plt.plot(w, P(w, *spirrid.get_samples(n)[:, i]),
                      color = 'black', label = 'random filament response')
        plt.plot(w, P(w, *spirrid.get_samples(n)[:, i]), color = 'black')
    plt.plot(w, spirrid.mu_q_arr, lw = 4, color = 'black', ls = '--',
              label = 'normalized yarn repsonse')
    plt.xlabel('$\mathrm{crack width} \, w \mathrm{[mm]}$', fontsize = 24)
    plt.ylabel('$\mathrm{force} \, P_\mathrm{y,0}/N_\mathrm{f} \mathrm{[N]}$', fontsize = 24)
    plt.xticks(fontsize = 20)
    plt.yticks(fontsize = 20)
    #plt.title( '$\mathrm{yarn \, crack \, bridge}$' , fontsize = 20 )
    plt.legend(loc = 'best')
    plt.show()

def random_samples_profiles(n):
    w = np.linspace(0, .5, 2)
    x = np.linspace(-50., 20., 100)
    P = CBEMClampedFiberSP()
    spirrid = SPIRRID(q = P,
                      sampling_type = 'LHS',
                      evars = dict(w = w,
                                    x = x),
                      tvars = dict(Ll = Ll,
                                   Lr = Lr,
                                   tau = RV('uniform', 0.0, .2),
                                   l = RV('uniform', 5.0, 15.0),
                                   A_f = Af,
                                   E_f = Ef,
                                   theta = RV('uniform', 0.0, .02),
                                   xi = xi,#RV('weibull_min', scale = 0.0179, shape = 5, n_int = 10),
                                   phi = phi,
                                   E_m = Em,
                                   A_m = Am,
                                   Nf = Nf
                                    ),
                    n_int = 20)

    for i in range(n):
        if i == n - 1:
            plt.plot(x, P(0.0, x, *spirrid.get_samples(n)[:, i]),
                      color = 'black', label = 'random filament response')
        plt.plot(x, P(0.0, x, *spirrid.get_samples(n)[:, i]), color = 'black')
    plt.plot(x, spirrid.mu_q_arr[1, :], lw = 4, color = 'black', ls = '--',
              label = 'normalized yarn repsonse')
    #plt.xlabel( '$\mathrm{position}{[mm]}$', fontsize = 24 )
    #plt.ylabel( '$\mathrm{force} \, P_\mathrm{y}/N_\mathrm{f} \mathrm{[N]}$', fontsize = 24 )
    #plt.title( '$\mathrm{yarn \, crack \, bridge}$' , fontsize = 20 )
    plt.xticks(fontsize = 20)
    plt.yticks(fontsize = 20)
    plt.legend(loc = 'lower right')
    plt.ylim(0, 0.35)
    plt.show()

def mlab_plot():
    w = np.linspace(0, 1.8, 30)
    x = np.linspace(-50., 20., 30)
    P = CBEMClampedFiberSP()
    spirrid = SPIRRID(q = P,
                      sampling_type = 'LHS',
                      evars = dict(w = w,
                                    x = x),
                      tvars = dict(Ll = 50.,
                                   Lr = 20.,
                                   tau = RV('uniform', 0.05, .15),
                                   l = RV('uniform', 2.0, 15.0),
                                   A_f = Af,
                                   E_f = Ef,
                                   theta = 0.01, #RV('uniform', 0.0, .02),
                                   xi = RV('weibull_min', scale = 0.0179, shape = 5, n_int = 10),
                                   phi = phi,
                                   E_m = Em,
                                   A_m = Am,
                                   Nf = 1.
                                    ),
                    n_int = 20)

    e_arr = make_ogrid([x, w])
    n_e_arr = [ e / np.max(np.fabs(e)) for e in e_arr ]
    mu_q_arr = spirrid.mu_q_arr
    n_mu_q_arr = mu_q_arr / np.max(np.fabs(mu_q_arr))
    m.surf(n_e_arr[0], n_e_arr[1], n_mu_q_arr)
    from numpy import array

# ------------------------------------------- 
    scene = engine.scenes[0]
    scene.scene.background = (1.0, 1.0, 1.0)
    scene.scene.camera.position = [1.5028781189276834, 3.5681173520859848, 1.9543753549631095]
    scene.scene.camera.focal_point = [-0.29999999701976776, 0.5, 0.5]
    scene.scene.camera.view_angle = 30.0
    scene.scene.camera.view_up = [-0.14835983161589669, -0.35107055575000001, 0.92452086252733578]
    scene.scene.camera.clipping_range = [1.9820957553309271, 6.1977861427731007]
    scene.scene.camera.compute_view_plane_normal()
    module_manager = engine.scenes[0].children[0].children[0].children[0].children[0]
    module_manager.scalar_lut_manager.lut_mode = 'Greys'
    m.show()

def ux(F):

    Kf = Af * Ef
    Km = Am * Em
    p = 2 * np.sqrt(Af / pi) * pi
    T = tau * p
    a = F * Km / (T * (Kf + Km))
    x = np.linspace(0, Ll + l + 10, 700)
    u1 = (F - T * a) * x / Kf
    u2 = T * x ** 2 / 2. / Kf - (T * Ll - F) * x / Kf + 0.5 * T * (Ll - a) ** 2 / Kf
    u3 = F * x / Kf - F * Km * (2 * Ll - a) / 2. / (Km + Kf) / Kf
    ux = u1 * H(Ll - a - x) + u2 * H(x - Ll + a) * H(Ll - x) + u3 * H(x - Ll)
    ux = ux * H(Ll + l - x)
    plt.plot(x, u1, ls = '--', color = 'black')
    plt.plot(x, u2, ls = '--', color = 'black')
    plt.plot(x, u3, ls = '--', color = 'black')
    plt.plot(x[0:600], ux[0:600], lw = 2, color = 'black', label = 'u(x)')

    P = 1.5 * F
    u2d = T * x ** 2 / 2. / Kf - (T * Ll - P) * x / Kf
    u3d = P * x / Kf - T * Ll ** 2 / 2. / Kf
    uxd = u2d * H(Ll - x) + u3d * H(x - Ll)
    #plt.plot( x, u2d )
    #plt.plot( x, u3d )
    plt.plot(x[0:600], uxd[0:600], color = 'black', lw = 2)

    plt.legend(loc = 'best')
    plt.ylim(-0.1)
    plt.xticks(fontsize = 18)
    plt.yticks(fontsize = 18)
    plt.show()

def ex(F):

    Kf = Af * Ef
    Km = Am * Em
    p = 2 * np.sqrt(Af / pi) * pi
    T = tau * p
    a = F * Km / (T * (Kf + Km))
    x = np.linspace(0, Ll + l + 10, 700)
    e1 = (F - T * a) / Kf
    e2 = T * x / Kf - (T * Ll - F) / Kf
    e3 = F / Kf
    ex = e1 * H(Ll - a - x) + e2 * H(x - Ll + a) * H(Ll - x) + e3 * H(x - Ll)
    ex = ex * H(Ll + l - x)
    plt.plot(x, len(x) * [e1], ls = '--', color = 'black')
    plt.plot(x, e2, ls = '--', color = 'black')
    plt.plot(x, len(x) * [e3], ls = '--', color = 'black')
    plt.plot(x[0:600], ex[0:600], lw = 2, color = 'black', label = '$\epsilon(x)$')

    P = F * 1.5
    e2d = T * x / Kf - (T * Ll - P) / Kf
    e3d = P / Kf
    exd = e2d * H(Ll - x) + e3d * H(x - Ll)

    plt.plot(x[0:600], exd[0:600], lw = 2, color = 'black')

    plt.legend(loc = 'best')
    plt.xticks(fontsize = 18)
    plt.yticks(fontsize = 18)
    plt.show()

def Puw():
    Kf = Af * Ef
    Km = Am * Em
    p = 2 * np.sqrt(Af / pi) * pi
    T = tau * p
    m = Km / (Kf + Km)
    u = np.linspace(0, 0.15, 200)
    w = np.linspace(0, 0.15, 200)
    cu = (Ll + l / 2.) * T - m * Ll * T
    Pu = (np.sqrt(cu ** 2 + 2 * T * m ** 2 * u * Kf) - cu) / m ** 2
    cw = (Ll + Lr + l) * T - m * (Ll + Lr) * T
    Pw = (np.sqrt(cw ** 2 + 4 * T * m ** 2 * w * Kf) - cw) / m ** 2 / 2.

    plt.plot(u, Pu, color = 'black', lw = 2, ls = '--', label = 'PO')
    plt.plot(w, Pw, label = 'CB', color = 'black', lw = 2)
    plt.legend(loc = 'best')
    plt.xticks(fontsize = 14)
    plt.yticks(fontsize = 14)
    plt.show()

def r_field():
    from stats.misc.random_field.random_field_1D import RandomField
    rf = RandomField(xgrid = np.linspace(0, 1000, 500))
    rf.distribution = 'Weibull'    
    rf.lacor = 3.0
    rf.shape = 10.
    rf.scale = 5.
    rf.loc = 0.0
    x = rf.xgrid
    y = rf.random_field
    plt.plot(x,y, color = 'grey', lw = 1, label = 'lacor = %.1f' %rf.lacor)
    rf.lacor = 10.
    y = rf.random_field
    plt.plot(x,y, color = 'black', lw = 2, label = 'lacor = %.1f' %rf.lacor)
    plt.legend(loc = 'best')
    plt.xticks(fontsize = 20)
    plt.yticks(fontsize = 20)
    plt.show()

def lacor():
    from stats.misc.random_field.random_field_1D import RandomField
    rf = RandomField()
    x = np.linspace(0,30, 200)
    y = rf.acor(x,10.)
    yy = rf.acor(x,3.0)
    plt.plot(x, y, color = 'black', lw = 2)
    plt.plot(x, yy, color = 'black', lw = 2)
    plt.xticks(fontsize = 20)
    plt.yticks(fontsize = 20)
    plt.show()

def cracks():
    from quaducom.meso.ctt.scm_numerical.ctt2 import CTT
    from stats.spirrid.spirrid import FunctionRandomization
    from stats.misc.random_field.random_field_1D import RandomField
    
    rf = CBEMClampedFiberSP()
    rand = FunctionRandomization(q = rf,
         evars = dict(w = np.linspace(0, .5, 40),
                       x = np.linspace(-20., 20., 101),
                       Ll = np.linspace(0.01, 20., 14),
                       Lr = np.linspace(0.01, 20., 14)),
         tvars = dict(tau = RV('uniform', loc = 0.5, scale = 0.2),
                      l = RV('uniform', loc = 0.0, scale = 2.0), A_f = 5.31e-4, E_f = 72e3, theta = 0.0,
                       xi = 1e20, phi = 1.0, E_m = 30e3, A_m = 50., Nf = 1700.),
                                n_int = 10)

    rand_field = RandomField( lacor = 1.0, xgrid = np.linspace(0, 60, 500),
                        nsim = 1, loc = 0.0, shape = 6.0, scale = 5.0,
                        non_negative_check = True, distribution = 'Weibull')
    
    ctt = CTT(length = 60., nx = 300,
          random_field = rand_field,
          cb_randomization = rand, cb_type = 'mean',
          force_min = 0.1, force_max = 400, n_force = 500)
    ctt.evaluate()
    plt.figure()
    plt.plot(ctt.x_arr, ctt.matrix_strength[-1,:], lw = 1, color = 'black')
    plt.plot(ctt.x_arr, ctt.sigma_m[200,:], lw = 2, color = 'black')
    plt.xticks(fontsize = 20)
    plt.yticks(fontsize = 20)
    plt.ylim(0,8)
    plt.show()
    plt.figure()
    plt.plot(ctt.x_arr, ctt.matrix_strength[-1,:], lw = 1, color = 'black')
    plt.plot(ctt.x_arr, ctt.sigma_m[250,:], lw = 2, color = 'black')
    plt.xticks(fontsize = 20)
    plt.yticks(fontsize = 20)
    plt.ylim(0,8)
    plt.show()
    plt.figure()
    plt.plot(ctt.x_arr, ctt.matrix_strength[-1,:], lw = 1, color = 'black')
    plt.plot(ctt.x_arr, ctt.sigma_m[450,:], lw = 2, color = 'black')
    plt.xticks(fontsize = 20)
    plt.yticks(fontsize = 20)
    plt.ylim(0,8)
    plt.show()  

# three phases of filament crack bridge
#Pw()
# force profile along a filament
#SP()
# "random" samples of a filament and average response l-d curves
random_samples_Pw(4)
# "random" samples of a filament and average response profiles
#random_samples_profiles(4)
# 3D plot
#mlab_plot()
# displacement along x in a fiber pullout at the force F ux(F)
#ux( .3 )
# strains along x in a fiber pullout at the force F ex(F)
#ex( .3 )
# PO vs CB
#Puw()
# random field
#r_field()
# autocorrelation function
#lacor()
# cracking profiles
#cracks()
