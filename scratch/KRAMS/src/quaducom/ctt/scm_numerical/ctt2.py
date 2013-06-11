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
# Created on Oct 21, 2011 by: rch

from enthought.traits.api import \
    HasTraits, Float, List, Range, \
    Callable, Instance, Trait, Int, \
    cached_property, Property, \
    implements, Interface, Array, WeakRef, \
    DelegatesTo, on_trait_change

from stats.spirrid.spirrid import SPIRRID, FunctionRandomization, MonteCarlo
from stats.spirrid.rv import RV
from quaducom.ctt.scm_numerical.interpolated_spirrid import InterpolatedSPIRRID
from stats.misc.random_field.gauss_1D import GaussRandomField
from operator import attrgetter

import numpy as np

class CBRandomSample(HasTraits):

    randomization = Instance(FunctionRandomization)

    sampling = Property
    @cached_property
    def _get_sampling(self):
        return MonteCarlo(randomization = self.randomization)

    def __call__(self, e):
        q = self.randomization.q
        sig = q(e, *self.sampling.theta)
        return np.sum(sig)

class CB(HasTraits):

    get_force_x_reinf = Callable
    randomization = WeakRef(FunctionRandomization)
    q = DelegatesTo('randomization')
    position = Float
    Ll = Float
    Lr = Float
    x = Array

    def get_eps_x_reinf(self, P):
        '''
        evaluation of strain profile in the vicinity of a crack bridge
        '''
        return self.get_force_x_reinf(P, self.x, self.Ll, self.Lr) / self.q.Kr

    def get_sigma_x_reinf(self, P):
        '''
        evaluation of stress profile in the vicinity of a crack bridge
        '''
        return self.get_force_x_reinf(P, self.x, self.Ll, self.Lr) / self.q.Ar

    def get_eps_x_matrix(self, P):
        '''
        evaluation of stress profile in the vicinity of a crack bridge
        '''
        Pff = np.ones(len(self.x)) * P
        return (Pff - self.get_force_x_reinf(P, self.x, self.Ll, self.Lr)) / self.q.Km

    def get_sigma_x_matrix(self, P):
        '''
        evaluation of stress profile in the vicinity of a crack bridge
        '''
        return self.get_eps_x_matrix(P) * self.q.E_m

class ICBMFactory(Interface):

    def new_cb_model(self):
        pass

class CBMFactory(HasTraits):

    implements(ICBMFactory)

    randomization = Instance(FunctionRandomization)

class CBMeanFactory(CBMFactory):

    # homogenization model
    spirrid = Property(Instance(SPIRRID))
    @cached_property
    def _get_spirrid(self):
        args = self.randomization.trait_get(['q', 'evars', 'tvars'])
        return SPIRRID(**args)

    interpolated_spirrid = Property()
    @cached_property
    def _get_interpolated_spirrid(self):
        return InterpolatedSPIRRID(spirrid = self.spirrid)

    #===========================================================================
    # Construct new crack bridge (the mean response is shared by all instances)
    #===========================================================================
    def new_cb(self):
        return CB(randomization = self.randomization,
                  get_force_x_reinf = self.interpolated_spirrid
                  )

class CBRandomFactory(CBMFactory):

    #===========================================================================
    # Construct new crack bridge (each crack bridge has its own random realization)
    #===========================================================================
    def new_cb(self):
        sample = CBRandomSample(randomization = self.randomization)
        return CB(get_force_x_reinf = sample,
                  randomization = self.randomization)

class CTT(HasTraits):

    cb_randomization = Instance(FunctionRandomization)

    cb_type = Trait('mean', dict(mean = CBMeanFactory,
                                 random = CBRandomFactory))

    cb_factory = Property(depends_on = 'cb_type')
    @cached_property
    def _get_cb_factory(self):
        return self.cb_type_(randomization = self.cb_randomization)

    cb_list = List

    length = Float(desc = 'composite specimen length')

    nx = Int(desc = 'number of discretization points')

    applied_force = Float

    x_arr = Property(Array, depends_on = 'length, nx')
    @cached_property
    def _get_x_arr(self):
        '''discretizes the specimen length'''
        return np.linspace(0, self.length, self.nx)

    sigma_mx_ff = Property(depends_on = 'applied_force, cb_randomization')
    @cached_property
    def _get_sigma_mx_ff(self):
        '''stress in the matrix along an uncracked composite - far field matrix stress'''
        sig_ff = self.applied_force * self.cb_randomization.q.Km / self.cb_randomization.q.Kc / self.cb_randomization.q.A_m
        return np.ones(len(self.x_arr)) * sig_ff

    matrix_strength = Array(desc = 'random matrix strength')

    def sort_cbs(self):
        '''sorts the CBs by position and adjusts the boundary conditions'''
        # sort the CBs
        self.cb_list = sorted(self.cb_list, key = attrgetter('position'))
        # specify the boundaries at the ends (0 and self.length) of the specimen 
        self.cb_list[0].Ll = self.cb_list[0].position
        self.cb_list[-1].Lr = self.length - self.cb_list[-1].position

        # specify the boundaries between the cracks
        for i, cb in enumerate(self.cb_list[:-1]):
            self.cb_list[i].Lr = (self.cb_list[i + 1].position - cb.position) / 2.
        for i, cb in enumerate(self.cb_list[1:]):
            self.cb_list[i + 1].Ll = (cb.position - self.cb_list[i].position) / 2.

        # specify the x range within the specimen length for every crack
        for i, cb in enumerate(self.cb_list):
            mask1 = self.x_arr > (cb.position - cb.Ll)
            if i == 0:
                mask1[0] = True
            mask2 = self.x_arr <= (cb.position + cb.Lr)
            cb.x = self.x_arr[mask1 * mask2] - cb.position

    @on_trait_change('applied_force')
    def insert_next_cb(self):
        '''seek for the minimum strength redundancy to find the position
        of the next crack'''

        while np.sum(self.sigma_m_x >= self.matrix_strength) > 0.:
            cr_pos = np.argmax(self.sigma_m_x - self.matrix_strength)
            cbf = self.cb_factory
            new_cb = cbf.new_cb()
            new_cb.position = self.x_arr[cr_pos]
            self.cb_list.append(new_cb)
            self.sort_cbs()
        # integrates the strain profile of the reinforcement and adds it to a list
        self.eps.append(np.trapz(self.eps_r_x, self.x_arr) / self.length)

    sigma_m_x = Property(depends_on = 'applied_force, cb_list')
    @cached_property
    def _get_sigma_m_x(self):
        if len(self.cb_list) == 0:
            profile = self.sigma_mx_ff
        else:
            profile = np.array([])
            for cb in self.cb_list:
                profile = np.hstack((profile, cb.get_sigma_x_matrix(self.applied_force)))
        #plt.plot(self.x_arr, profile)
        return profile

    eps_m_x = Property(depends_on = 'applied_force, cb_list')
    @cached_property
    def _get_eps_m_x(self):
        return self.sigma_m_x / self.cb_randomization.q.E_m

    eps_r_x = Property(depends_on = 'applied_force, cb_list')
    @cached_property
    def _get_eps_r_x(self):
        if len(self.cb_list) == 0:
            profile = self.sigma_mx_ff / self.cb_randomization.q.E_m
        else:
            profile = np.array([])
            for cb in self.cb_list:
                profile = np.hstack((profile, cb.get_eps_x_reinf(self.applied_force)))
        return profile

    eps = List

    sig = Property
    @cached_property
    def _get_sig(self):

        cbf = self.cb_factory
        total = 0.0

        for i in range(0, 10):
            new_cb = cbf.new_cb()
            new_cb.set(position = i)
            self.cb_list.append(new_cb)
            # get the response
            for cb in self.cb_list:
                total += cb.get_force_x_reinf(i)

        return total

if __name__ == '__main__':

    from quaducom.resp_func.cb_emtrx_clamped_fiber import \
        CBEMClampedFiberSP
    import enthought.mayavi.mlab as m
    from stats.spirrid import orthogonalize
    from matplotlib import pyplot as plt

    length = 1000
    nx = 200
    matrix_strength = GaussRandomField(lacor = 5.,
                                    xgrid = np.linspace(0, length, nx),
                                    nsim = 1,
                                    mean = 4.,
                                    stdev = 1.,
                                    non_negative_check = True
                               ).random_field

    rf = CBEMClampedFiberSP()
    rand = FunctionRandomization(q = rf,
         evars = dict(w = np.linspace(0.0, 1.0, 20),
                       x = np.linspace(-60., 60., 30),
                       Ll = np.linspace(1., 100., 5),
                       Lr = np.linspace(1., 100., 5),
                        ),
         tvars = dict(tau = 0.15, #RV( 'uniform', 0.7, 1.0 ),
                       l = RV('uniform', 2.0, 10.0),
                       A_r = 0.89,
                       E_f = 72e3,
                       theta = RV('uniform', 0.0, .10),
                       xi = 0.5, #RV( 'weibull_min', scale = 0.017, shape = 5, n_int = 10 ),
                       phi = 1.0,
                       E_m = 30e3,
                       A_m = rf.A_m,
                       Nf = 1700.
                        ),
         n_int = 20)

    ctt = CTT(length = length,
              nx = nx,
              matrix_strength = matrix_strength,
              cb_randomization = rand,
              cb_type = 'mean',
              )
    P = np.linspace(0.1, 400, 200)
    matrix = []
    reinf = []
    for p in P:
        ctt.applied_force = p
        matrix.append(ctt.eps_r_x)
        reinf.append(ctt.eps_m_x)
    e_arr = orthogonalize([ctt.x_arr, P])
    n_e_arr = [ e / np.max(np.fabs(e)) for e in e_arr ]

    mu_m_arr = np.array(matrix)
    mu_r_arr = np.array(reinf)

    n_mu_m_arr = mu_m_arr / np.max(np.fabs(mu_m_arr))
    n_mu_r_arr = mu_r_arr / np.max(np.fabs(mu_r_arr))

    #m.surf(n_e_arr[0], n_e_arr[1], n_mu_m_arr)
    m.surf(n_e_arr[0], n_e_arr[1], n_mu_r_arr)
    plt.plot(ctt.eps, P)
    plt.show()
    m.show()
