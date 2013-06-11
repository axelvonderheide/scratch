'''
Created on Jul 26, 2012

@author: rostar
'''

from etsproxy.traits.api import \
    Instance, Array, List, cached_property, Property, Float
from etsproxy.traits.ui.api import ModelView
from stats.spirrid.rv import RV
from stats.misc.random_field.random_field_1D import RandomField
import numpy as np
from matplotlib import pyplot as plt
from scratch.scm.scm_dependent import SCM


class SCMView(ModelView):

    model = Instance(SCM)

    def crack_widths(self, sigma_c):
        # find the index of the nearest value in the load range
        idx = np.abs(self.model.load_sigma_c - sigma_c).argmin()
        # evaluate the relative strain e_rel between fibers
        # and matrix for the given load
        e_rel = self.eps_f[idx, :] - self.eps_m[idx, :]
        # pick the cracks that emerged at the given load
        cb_load = self.model.cb_list(sigma_c)
        if cb_load[0] is not None:
            # find the symmetry points between cracks as
            # the 0 element of their x range
            indexes = []
            for cb in cb_load:
                indexes.append(np.where(cb.position +
                                cb.x[0] == self.model.x_arr)[0])
            # add the index of the last point
            indexes.append(self.model.nx - 1)
            # list of crack widths to be filled in a loop with integrated e_rel
            crack_widths = [np.trapz(e_rel[idx:indexes[i + 1]],
                            self.model.x_arr[idx:indexes[i + 1]])
                            for i, idx in enumerate(indexes[:-1])]
            return np.array(crack_widths)
        else:
            return 0.0

    eval_w = Property(List, depends_on='model')

    @cached_property
    def _get_eval_w(self):
        return [self.crack_widths(load) for load in self.model.load_sigma_c]

    w_mean = Property(Array, depends_on='model')

    @cached_property
    def _get_w_mean(self):
        return np.array([np.mean(w) for w in self.eval_w])

    w_median = Property(Array, depends_on='model')

    @cached_property
    def _get_w_median(self):
        return np.array([np.median(w) for w in self.eval_w])

    w_stdev = Property(Array, depends_on='model')

    @cached_property
    def _get_w_stdev(self):
        return np.array([np.std(w) for w in self.eval_w])

    w_max = Property(Array, depends_on='model')

    @cached_property
    def _get_w_max(self):
        return np.array([np.max(w) for w in self.eval_w])

    x_area = Property(depends_on='model.')

    @cached_property
    def _get_x_area(self):
        return  np.ones_like(self.model.load_sigma_c)[:, np.newaxis] \
            * self.model.x_arr[np.newaxis, :]

    eps_m = Property(Array, depends_on='model.')

    @cached_property
    def _get_eps_m(self):
        return self.model.sigma_m_x / self.model.cb_randomization.tvars['E_m']

    sigma_f = Property(Array, depends_on='model.')
    @cached_property
    def _get_sigma_f(self):
        V_f = self.model.cb_randomization.tvars['V_f']
        return (self.model.load_sigma_c[:, np.newaxis] -
                self.model.sigma_m_x * (1. - V_f)) / V_f

    eps_f = Property(Array, depends_on='model.')
    @cached_property
    def _get_eps_f(self):
        return self.sigma_f / self.model.cb_randomization.tvars['E_f']

    eps_sigma = Property(depends_on='model.')
    @cached_property
    def _get_eps_sigma(self):
        eps = np.trapz(self.eps_f, self.x_area, axis=1) / self.model.length
        eps = eps[np.isnan(eps) == False]
        if len(eps) != len(self.model.load_sigma_c):
            eps = list(eps) + [list(eps)[-1]]
            sigma = self.model.load_sigma_c[:len(eps)]
            sigma[-1] = 0.0
            return eps, sigma
        else:
            return eps, self.model.load_sigma_c

    def get_epsc(self, load):
        if len(self.model.cracks_list) is not 0:
            idx = np.sum(np.array(self.model.sigma_c_crack) < load) - 1
            if idx == -1:
                cb_list = [None]
            else:
                print self.model.cracks_list
                cb_list = self.model.cracks_list[idx]
        else:
            cb_list = [None]
        if cb_list == [None]:
            epsc = load / self.model.E_c
        else:
            u = 0.0
            for cb in cb_list:
                cb.load = load
                u += np.trapz(cb.get_epsf_x(), cb.x)
                plt.plot(cb.x, cb.get_epsf_x())
                plt.plot([cb.x[0], cb.x[-1]], [load / self.model.E_c,load / self.model.E_c])
                plt.show()
            epsc = u / self.model.length
        return epsc

if __name__ == '__main__':
    from scratch.CB.dependent_fibers.reinforcement import Reinforcement
    from stats.spirrid.rv import RV
    from scratch.CB.dependent_fibers.composite_CB_model import CompositeCrackBridge
    from scratch.CB.dependent_fibers.reinforcement import WeibullFibers

    reinf = Reinforcement(r=0.00345,#RV('uniform', loc=0.002, scale=0.002),
                          tau=RV('uniform', loc=.1, scale=.1),
                          V_f=0.1,
                          E_f=200e3,
                          xi=100.03,
                          n_int=20)
    E_m = 25e3
    length = 20.
    nx = 500
    random_field = RandomField(seed=True,
                               lacor=1.,
                               xgrid=np.linspace(0., length, 200),
                               nsim=1,
                               loc=.0,
                               shape=15.,
                               scale=6.,
                               non_negative_check=True,
                               distribution='Weibull'
                               )
    scm = SCM(length=length,
              nx=nx,
              random_field=random_field,
              E_m=E_m,
              reinforcement_lst=[reinf],
              load_sigma_c_min=.1,
              load_sigma_c_max=11.,
              load_n_sigma_c=200,
              )
    scm_view = SCMView(model=scm)
    scm_view.model.evaluate()

    def plot():
        load = np.linspace(1., 12., 30)
        epsc = []
        for l in load:
            epsc.append(scm_view.get_epsc(l))
        plt.plot(epsc, load, color='black', lw=2, label='model')
        plt.legend(loc='best')
        plt.xlabel('eps_c [mm]')
        plt.ylabel('composite stress [MPa]')
#        plt.figure()
#        plt.hist(scm_view.crack_widths(16.), bins=20, label='load = 20 MPa')
#        plt.hist(scm_view.crack_widths(13.), bins=20, label='load = 15 MPa')
#        plt.hist(scm_view.crack_widths(10.), bins=20, label='load = 10 MPa')
#        plt.legend(loc='best')
#        plt.figure()
#        plt.plot(scm.load_sigma_c, scm_view.w_mean,
#                 color='green', lw=2, label='mean crack width')
#        plt.plot(scm.load_sigma_c, scm_view.w_median,
#                 color='blue', lw=2, label='median crack width')
#        plt.plot(scm.load_sigma_c, scm_view.w_mean + scm_view.w_stdev,
#                 color='black', label='stdev')
#        plt.plot(scm.load_sigma_c, scm_view.w_mean - scm_view.w_stdev,
#                 color='black')
#        plt.plot(scm.load_sigma_c, scm_view.w_max,
#                 ls='dashed', color='red', label='max crack width')
#        plt.legend(loc='best')
        plt.show()
    plot()
