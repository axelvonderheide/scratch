from enthought.traits.api import \
    Float, Str, implements

from enthought.traits.ui.ui_traits import Image

from enthought.traits.ui.menu import OKButton, CancelButton

from enthought.traits.ui.api import \
    View, Item
from quaducom.resp_func.cb_short_fiber import CBShortFiber
from math import e, pi
from numpy import sqrt, linspace, sign, abs, cos
from stats.spirrid.i_rf import IRF
from stats.spirrid.rf import RF

from matplotlib import pyplot as plt

def H(x):
    return sign(sign(x) + 1.)

class TESTCBShortFiber(RF):
    '''
    Crack bridged by a short fiber with constant
    frictional interface to the linear elastic matrix''
    '''

    implements(IRF)

    E_f = Float(200e+3 , auto_set=False, enter_set=True,
                desc='filament stiffness [N/mm2]',
                distr=['uniform', 'norm'],
                scale=210e3, shape=0)

    D_f = Float(0.3, auto_set=False, enter_set=True,
                desc='filament diameter [mm]',
                distr=['uniform', 'norm'],
                scale=0.5, shape=0)

    L_e = Float(8.5, auto_set=False, enter_set=True,
                desc='shorter embedded length [mm]',
                distr=['uniform'],
                scale=8.5, shape=0)

    L_f = Float(17.0, auto_set=False, enter_set=True,
                desc='fiber length [mm]',
                distr=['uniform', 'norm'],
                scale=30, shape=0)

    tau = Float(1.76, auto_set=False, enter_set=True,
                desc='bond shear stress [N/mm2]',
                distr=['norm', 'uniform'],
                scale=1.76, shape=0.5)

    f = Float(0.03, auto_set=False, enter_set=True,
            desc='ging coefficient',
            distr=['uniform', 'norm'],
                scale=0.05, shape=0)
    

    phi = Float(0.0, auto_set=False, enter_set=True,
       desc='inclination angle',
       distr=['sin2x', 'sin_distr'],
                scale=1.0, shape=0)

    l = Float(0.0, auto_set=False, enter_set=True,
              distr=['uniform'], desc='free length')

    nu = Float(1.0 , desc='matrix/composite ratio' , auto_set=False, enter_set=True,
                distr=['uniform', 'norm'],
                scale=0.5, shape=0)
    
    w = Float(ctrl_range=(0, 0.01, 100), auto_set=False, enter_set=True)

    x_label = Str('crack_opening [mm]', enter_set=True, auto_set=False)
    y_label = Str('force [N]', enter_set=True, auto_set=False)


    def __call__(self, w, tau, L_f, D_f, E_f, L_e, phi, f , nu, L_cr):

        w = w 
        T = tau * pi * D_f
        A_f = D_f ** 2 / 4. * pi
        cA = L_cr * T * (1 - nu)
        g = e ** (f * phi)
        # debonding stage
        q_deb = 0.5 / nu ** 2 * (sqrt(cA ** 2 + 4 * w * E_f * A_f * nu ** 2 * T) - cA) * g
        
        # displacement at which debonding is finished
        w0 = L_e * (cA + g * L_e * T * nu ** 2) / (E_f * A_f)
        #q_w0=0.5 / nu ** 2 * (sqrt(cA ** 2 + 4 * w * E_f * A_f * nu ** 2 * T) - cA) * g
        # pulling out stage - the fiber is pulled out from the
        # side with the shorter embedded length only
        q_pull = L_e * T * g * ((w0 - w) / ((L_e + 1e-15) - w0) + 1)
        q = q_deb * H(L_e * T * g - q_deb) + q_pull * H(q_deb - L_e * T * g)

        # include inclination influence
        q = q * H(q) 

        
        return q


if __name__ == '__main__':
    import numpy as np
    ############ela##########
    q = TESTCBShortFiber()
    w = np.linspace(0, 0.01, 20000)
    #          w,tau,L_f,D_f, E_f  , L_e, phi, f , nu , Lcr
    result_e = q(w, 2, 10, 0.3, 200000, 5, 1, 1 , 0.8, 20)
    ###########rigid#############
    s = CBShortFiber()
    #          w, tau, L_f, D_f, E_f, le, phi, f, l, theta, xi
    result_r = s(w, 2, 10, 0.3, 200000, 5, 1, 1 , 0, 0, 1e15)
    #result_r2 = s(w, 4, 10, 0.3, 200000, 5, 1, 1 , 0, 0, 1e15)
    #########Plot############
    plt.plot(w, result_r, 'black')
    #plt.plot(w, result_r2, 'black')
    plt.plot(w, result_e, 'red')
    plt.plot(w, abs(result_e - result_r), 'green')
    plt.show()
