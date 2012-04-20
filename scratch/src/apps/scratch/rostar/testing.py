import numpy as np
from math import e
from matplotlib import pyplot as plt
from scipy.stats import norm
from scipy.optimize import fsolve as op

from quaducom.resp_func.cb_clamped_fiber import H
from quaducom.resp_func.cb_clamped_fiber import CBClampedFiber, CBClampedFiberSP
from quaducom.resp_func.cb_emtrx_clamped_fiber import CBEMClampedFiber, CBEMClampedFiberSP

w = np.linspace( 0, .2, 300 )
x = np.linspace( -60, 30, 300 )
def Pw():

    Pe = CBEMClampedFiber()
    qe = Pe( w, 6.37, 10., 2., 72e3, 30000., 50., 0.0, 999, 1., 15., 30. )
    P = CBClampedFiber()
    q = P( w, 6.37, 10., 2., 72e3, 0.0, 999, 1., 15., 30. )
    plt.plot( w, q, label = 'rigid', color = 'black', lw = 2 )
    plt.plot( w, qe, label = 'elastic', color = 'blue', lw = 2 )
    plt.legend( loc = 'best' )
    plt.show()

def SP():
    spe = CBEMClampedFiberSP()
    qe = spe( .02, x, 2.37, 10., 3., 72e3, 30000., 50., 0.0, 999, 1., 50., 30. )
    sp = CBClampedFiberSP()
    q = sp( .02, x, 2.37, 10., 3., 72e3, 0.0, 999, 1., 50., 30. )
    plt.plot( x, q, lw = 2, color = 'black', label = 'rigid' )
    plt.plot( x, qe, lw = 2, color = 'blue', label = 'elastic' )
    plt.xticks( fontsize = 14 )
    plt.yticks( fontsize = 14 )
    plt.legend( loc = 'best' )
    plt.show()

Pw()
SP()
