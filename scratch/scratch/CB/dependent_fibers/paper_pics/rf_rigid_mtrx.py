'''
Created on Nov 20, 2012

@author: rostar
'''
#-------------------------------------------------------------------------------
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

from etsproxy.traits.api import \
    Float, Str, implements, cached_property, Property
import numpy as np
from stats.spirrid.i_rf import \
    IRF
from stats.spirrid.rf import \
    RF
from matplotlib import pyplot as plt


def H(x):
    return x >= 0.0


class CBRigidMatrix(RF):

    implements(IRF)

    title = Str('crack bridge - clamped fiber with constant friction')

    xi = Float(0.0179, auto_set=False, enter_set=True, input=True,
                distr=['weibull_min', 'uniform'])

    tau = Float(0.5, auto_set=False, enter_set=True, input=True,
                distr=['uniform', 'norm'])

    E_f = Float(72.0e3, auto_set=False, enter_set=True, input=True,
                  distr=['uniform'])

    V_f = Float(0.2, auto_set=False, enter_set=True, input=True,
                  distr=['uniform'])

    r = Float(0.00345, auto_set=False, enter_set=True, input=True,
                  distr=['norm', 'uniform'])

    w = Float(auto_set=False, enter_set=True, input=True,
               distr=['uniform'], desc='crack width',
               ctrl_range=(0.0, 1.0, 10))

    C_code = Str('')

    def __call__(self, w, r, tau, E_f, V_f, xi):
        T = 2. * tau / r
        eps = np.sqrt(T * w / E_f)
        # include breaking strain
        eps = eps * H(xi - eps)
        return eps


class CBRigidMatrixSP(CBRigidMatrix):
    x = Float(0.0, auto_set=False, enter_set=True, input=True,
              distr=['uniform'], desc='distance from crack')
    C_code = Str('')

    def __call__(self, w, x, r, tau, E_f, V_f, xi):
        T = 2. * tau / r
        epsf0 = super(CBRigidMatrixSP, self).__call__(w, r, tau, E_f, V_f, xi)
        epsfx = epsf0 - (T / E_f * abs(x))
        eps = np.sqrt(T * w / E_f)
        survived_fibers = np.sum(H(xi - eps))/float(len(epsf0.flatten()))
        print survived_fibers
        return epsfx * H(epsfx) / survived_fibers

if __name__ == '__main__':
    
    def sigma_w():
        cb = CBRigidMatrix()
        w = np.linspace(0.,0.5, 100)
        plt.plot(w, 200e3*cb(w, 0.00345, 0.5, 200e3, 0.3, 100.))
        plt.show()
        
    def SP():
        cb = CBRigidMatrixSP()
        x = np.linspace(0., 2., 100)
        plt.plot(x, cb(0.01, x, 0.002, 0.5, 200e3, 0.3, 100.))
        plt.show()
    SP()
