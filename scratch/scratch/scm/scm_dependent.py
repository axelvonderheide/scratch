
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
# Created on Oct 21, 2011 by: rch

from etsproxy.traits.api import \
    HasTraits, Instance, Int, Array, List, \
    implements, Trait, cached_property, Property, Float
from scratch.scm.interpolated_response import InterpolatedCB
from stats.misc.random_field.random_field_1D import RandomField
from operator import attrgetter
import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import brentq
import copy
from scratch.CB.dependent_fibers.composite_CB_postprocessor import \
    CompositeCrackBridgeView
from matplotlib import pyplot as plt
from scratch.CB.dependent_fibers.reinforcement import Reinforcement
from scratch.CB.dependent_fibers.composite_CB_model import CompositeCrackBridge


class SCM(HasTraits):
    '''Stochastic Cracking Model - compares matrix strength and stress,
    inserts new CS instances at positions, where the matrix strength
    is lower than the stress; evaluates stress-strain diagram
    by integrating the strain profile along the composite'''

    E_m = Float
    reinforcement_lst = List(Instance(Reinforcement))
    length = Float(desc='composite specimen length')
    nx = Int(desc='number of discretization points')
    sigma_c_crack = List
    cracks_list = List

    load_sigma_c = Property(depends_on='+load')
    @cached_property
    def _get_load_sigma_c(self):
        #applied external load in terms of composite stress
        return np.linspace(self.load_sigma_c_min,
                           self.load_sigma_c_max, self.load_n_sigma_c)

    load_sigma_c_min = Float(load=True)
    load_sigma_c_max = Float(load=True)
    load_n_sigma_c = Int(load=True)

    x_arr = Property(Array, depends_on='length, nx')
    @cached_property
    def _get_x_arr(self):
        # discretizes the specimen length
        return np.linspace(0., self.length, self.nx)

    random_field = Instance(RandomField)
    matrix_strength = Property(depends_on='random_field.+modified')

    @cached_property
    def _get_matrix_strength(self):
        # evaluates a random field
        # realization and creates a spline reprezentation
        rf = self.random_field.random_field
        rf_spline = interp1d(self.random_field.xgrid, rf)
        return rf_spline(self.x_arr)

    def sort_cbs(self):
        # sorts the CBs by position and adjusts the boundary conditions
        # sort the CBs
        cb_list = self.cracks_list[-1]
        crack_position = cb_list[-1].position
        cb_list = sorted(cb_list, key=attrgetter('position'))
        #find idx of the new crack
        for i, crack in enumerate(cb_list):
            if crack.position == crack_position:
                idx = i
        # specify the boundaries
        if idx != 0:
            # there is a crack at the left hand side
            cbl = cb_list[idx - 1]
            cb = cb_list[idx]
            cbl.Lr = (cb.position - cbl.position) / 2.
            cb.Ll = cbl.Lr
        else:
            # the new crack is the first from the left hand side
            cb_list[idx].Ll = cb_list[idx].position

        if idx != len(cb_list) - 1:
            # there is a crack at the right hand side
            cb, cbr = cb_list[idx], cb_list[idx + 1]
            cbr.Ll = (cbr.position - cb.position) / 2.
            cb.Lr = cbr.Ll
        else:
            # the new crack is the first from the right hand side
            cb_list[idx].Lr = self.length - cb_list[idx].position

        # specify the x range and stress profile for
        # the new crack and its neighbors
        idxs = [idx - 1, idx, idx + 1]
        if idx == 0:
            idxs.remove(-1)
        if idx == len(cb_list) - 1:
            idxs.remove(len(cb_list))
        for idx in idxs:
            mask1 = self.x_arr >= (cb_list[idx].position - cb_list[idx].Ll)
            if idx == 0:
                mask1[0] = True
            mask2 = self.x_arr <= (cb_list[idx].position + cb_list[idx].Lr)
            cb_list[idx].x = self.x_arr[mask1 * mask2] - cb_list[idx].position
        self.cracks_list[-1] = cb_list

    def cb_list(self, load):
        if len(self.cracks_list) is not 0:
            idx = np.sum(np.array(self.sigma_c_crack) < load) - 1
            return self.cracks_list[idx]
        else:
            return [None]

    E_c = Property(depends_on='+reinforcement_list, E_m')
    @cached_property
    def _get_E_c(self):
        Kf = 0.0
        Vf_tot = 0.0
        for reinf in self.reinforcement_lst:
            Kf += reinf.V_f * reinf.E_f
            Vf_tot += reinf.V_f
        return self.E_m * (1. - Vf_tot) + Kf

    def sigma_m(self, load):
        sigma_m = load * self.E_m / self.E_c * np.ones(len(self.x_arr))
        cb_load = self.cb_list(load)
        if cb_load[0] is not None:
            for cb in cb_load:
                cb.load = load
                crack_position_idx = np.argwhere(self.x_arr == cb.position)
                left = crack_position_idx - len(np.nonzero(cb.x < 0.)[0])
                right = crack_position_idx + len(np.nonzero(cb.x > 0.)[0]) + 1
                sigma_m[left:right] = cb.get_sigmam_x().T
        return sigma_m

    def residuum(self, q):
        res = np.min(self.matrix_strength - self.sigma_m(q))
        return res

    def evaluate(self):
        # seek for the minimum strength redundancy to find the position
        # of the next crack
        last_pos = 0.0
        q_min = 0.0
        q_max = self.load_sigma_c_max
        while np.any(self.sigma_m(q_max) > self.matrix_strength):
            q_min = brentq(self.residuum, q_min, q_max)
            crack_position = self.x_arr[np.argmin(self.matrix_strength -
                                                  self.sigma_m(q_min))]
            model = CompositeCrackBridge(E_m=self.E_m,
                                         reinforcement_lst=self.reinforcement_lst)
            ccb_view = CompositeCrackBridgeView(model=model)
            new_cb = InterpolatedCB(cb=ccb_view)
            new_cb.position = float(crack_position)
            if len(self.cracks_list) is not 0:
                self.cracks_list.append(copy.copy(self.cracks_list[-1])
                                        + [new_cb])
            else:
                self.cracks_list.append([new_cb])
            self.sort_cbs()
            self.sigma_c_crack.append(q_min - self.load_sigma_c_max / 1e5)
#            plt.plot(self.x_arr, self.sigma_m(q_min))
#            plt.plot(self.x_arr, self.matrix_strength)
#            plt.show()
#            cb_list = self.cracks_list[-1]
#            cb = [CB for CB in cb_list if
#                  CB.position == float(crack_position)][0]
#            mu_q = cb.get_sigma_f_x_reinf(self.load_sigma_c,
#                                          np.array([0.0]),
#                                          cb.Ll, cb.Lr).flatten()
#            mu_q_real = mu_q[np.isnan(mu_q) == False]
#            new_q_max = np.max(mu_q_real) * self.cb_randomization.tvars['V_f']
#            new_q_max = self.load_sigma_c_max
#            if new_q_max < q_max:
#                q_max = new_q_max
            if float(crack_position) == last_pos:
                raise ValueError('''got stuck in loop,
                try to adapt x, w, BC ranges''')
            last_pos = float(crack_position)

    sigma_m_x = Property(depends_on='''random_field.+modified,
                            +load, nx, length, cb_type''')
    @cached_property
    def _get_sigma_m_x(self):
        sigma_m_x = np.zeros_like(self.load_sigma_c[:, np.newaxis]
                                  * self.x_arr[np.newaxis, :])
        for i, q in enumerate(self.load_sigma_c):
            sigma_m_x[i, :] = self.sigma_m(q)
        return sigma_m_x

if __name__ == '__main__':
    from stats.spirrid.rv import RV
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
              load_sigma_c_max=18.,
              load_n_sigma_c=200,
              )
    scm.evaluate()

