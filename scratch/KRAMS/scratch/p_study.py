'''
Created on Mar 6, 2012

@author: rostar
'''

from quaducom.meso.ctt.scm_numerical.scm_model import SCM
import etsproxy.mayavi.mlab as m
import numpy as np
from stats.spirrid import make_ogrid as orthogonalize
from matplotlib import pyplot as plt
from stats.spirrid.spirrid import FunctionRandomization
from stats.misc.random_field.random_field_1D import RandomField
from stats.spirrid.rv import RV

from quaducom.micro.resp_func.cb_emtrx_clamped_fiber_stress import \
    CBEMClampedFiberStressSP

if __name__ == '__main__':
    # filaments
    r = 0.00345
    Vf = 0.0103
    tau = .5 #RV('uniform', loc = 0.02, scale = .01) # 0.5
    Ef = 200e3
    Em = 25e3
    l = 0.0#RV( 'uniform', scale = 10., loc = 2. )
    theta = 0.0
    xi = RV( 'weibull_min', scale = 0.01, shape = 5 ) # 0.017
    phi = 1.

    length = 600.
    nx = 5000
    random_field = RandomField(seed = False,
                               lacor = 4.,
                                xgrid = np.linspace(0., length, 600),
                                nsim = 1,
                                loc = .0,
                                shape = 7.5,
                                scale = 7.8,
                                non_negative_check = True,
                                distribution = 'Weibull'
                               )

    rf = CBEMClampedFiberStressSP()
    rand = FunctionRandomization(   q = rf,
                                    tvars = dict(tau = tau,
                                                 l = l,
                                                 E_f = Ef,
                                                 theta = theta,
                                                 xi = xi,
                                                 phi = phi,
                                                 E_m = Em,
                                                 r = r,
                                                 V_f = Vf
                                                 ),
                                    n_int = 30
                                    )

    scm = SCM(length = length,
              nx = nx,
              random_field = random_field,
              cb_randomization = rand,
              cb_type = 'mean',
              load_sigma_c_min = 0.1,
              load_sigma_c_max = 30.,
              load_n_sigma_c = 100,
              n_w = 40,
              n_x = 101,
              n_BC = 4
              )
    
    label = ['discrete', 'random (COV = 23%%)', 'random (COV = 52%%)']
    def plot(i):
        eps, sigma = scm.eps_sigma
        mask = np.isnan(eps) == False
        eps = eps[mask]
        sigma = sigma[mask]
        sigma[-1] = 0.0
        plt.plot(eps, sigma, color = 'black', lw = 2, label = label[i])
        plt.xlabel('composite strain [-]')
        plt.ylabel('composite stress [MPa]')
    
    tau_lst = [0.1, .7]
                #RV('uniform', loc = 0.3, scale = .4),
                #RV('uniform', loc = 0.1, scale = .8)]
    for i,t in enumerate(tau_lst):
        scm = SCM(length = length,
              nx = nx,
              random_field = random_field,
              cb_randomization = rand,
              cb_type = 'mean',
              load_sigma_c_min = 0.1,
              load_sigma_c_max = 12.,
              load_n_sigma_c = 200,
              n_w = 50,
              n_x = 71,
              n_BC = 20
              )
        scm.cb_randomization.tvars['tau'] = t
        scm.evaluate()
        plot(i)
    
    plt.show()



