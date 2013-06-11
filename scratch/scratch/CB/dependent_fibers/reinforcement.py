
'''
Created on 23.10.2012

An instance or a list of instances of the Reinforcement class
can be used by the composite crack bridge model.

@author: Q
'''

import numpy as np
from stats.spirrid.rv import RV
from etsproxy.traits.api import HasTraits, cached_property, \
    Float, Property, Int, Str
from types import FloatType
from util.traits.either_type import EitherType
from scipy.stats import weibull_min
from math import pi

class WeibullFibers(HasTraits):
    '''class evaluating damage for weibull fibers with linearly decreasing stress'''
    shape = Float(5.0)
    scale = Float(0.02)
    L0 = Float(10.)
    
    distribution = Property(depends_on='shape, scale, L0')
    @cached_property
    def _get_distribution(self):
        return weibull_min(self.shape, scale=self.scale)
    
    def weibull_fibers_Pf(self, epsy_arr, depsf, x_short, x_long):
        x_short = np.hstack((x_short[1:], np.repeat(x_short[-1], len(epsy_arr)-len(x_short[1:]))))
        x_long = np.hstack((x_long[1:], np.repeat(x_long[-1], len(epsy_arr)-len(x_long[1:]))))
        Pf_short = (((depsf * x_short - 1.) * ((epsy_arr * (1. - depsf * x_short)))
                     / self.scale) ** self.shape + ( epsy_arr / self.scale ) ** self.shape) / (self.shape + 1) / depsf
        Pf_long = (((depsf * x_long - 1.) * ((epsy_arr * (1. - depsf * x_long)))
                     / self.scale) ** self.shape + ( epsy_arr / self.scale ) ** self.shape) / (self.shape + 1) / depsf
        Pf_short = 1. - np.exp(-1. / self.L0 * Pf_short)
        Pf_long = 1. - np.exp(-1. / self.L0 * Pf_long)
        return Pf_long + Pf_short - Pf_long * Pf_short  


class Reinforcement(HasTraits):

    label = Str('reinforcement')
    r = EitherType(klasses=[FloatType, RV])
    V_f = Float
    E_f = Float
    xi = EitherType(klasses=[FloatType, RV, WeibullFibers])
    tau = EitherType(klasses=[FloatType, RV])
    n_int = Int

    results = Property(depends_on='r, V_f, E_f, xi, tau, n_int')
    @cached_property
    def _get_results(self):
        stat_weights = 1.0
        if isinstance(self.tau, RV):
            tau = self.tau.ppf(
                np.linspace(.5 / self.n_int, 1. - .5 / self.n_int, self.n_int))
            stat_weights *= 1. / self.n_int
            nu_r_tau = np.ones_like(tau)
        else:
            tau = self.tau
            nu_r_tau = 1.0
        if isinstance(self.r, RV):
            r = self.r.ppf(
                np.linspace(.5 / self.n_int, 1. - .5 / self.n_int, self.n_int))
            stat_weights *= 1. / self.n_int
            r2 = r ** 2
            nu_r = r2 / np.mean(r2)
        else:
            r = self.r
            nu_r = nu_r_tau * 1.0
        if isinstance(tau, np.ndarray) and isinstance(r, np.ndarray):
            r = r.reshape(1, self.n_int)
            tau = tau.reshape(self.n_int, 1)
            nu_r_r = (r2 / np.mean(r2)).reshape(1, self.n_int)
            nu_r_tau = np.ones(self.n_int).reshape(self.n_int, 1)
            nu_r = nu_r_r * nu_r_tau
            return (2. * tau / r / self.E_f).flatten(), stat_weights, nu_r.flatten()
        else:
            return 2. * tau / r / self.E_f, stat_weights, nu_r

    depsf_arr = Property(depends_on='r, V_f, E_f, xi, tau, n_int')
    @cached_property
    def _get_depsf_arr(self):
        return self.results[0]

    stat_weights = Property(depends_on='r, V_f, E_f, xi, tau, n_int')
    @cached_property
    def _get_stat_weights(self):
        return self.results[1]
    
    nu_r = Property(depends_on='r, V_f, E_f, xi, tau, n_int')
    @cached_property
    def _get_nu_r(self):
        return self.results[2]  
    

