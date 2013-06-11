
'''
Created on 14.6.2012

@author: Q
'''

from stats.spirrid.sampling import FunctionRandomization
from stats.spirrid.rv import RV
from stats.misc.random_field.random_field_1D import RandomField
import numpy as np
from matplotlib import pyplot as plt
from quaducom.meso.ctt.scm_numerical.scm_model import SCM
from quaducom.micro.resp_func.cb_emtrx_clamped_fiber_stress import \
CBEMClampedFiberStressSP
from quaducom.meso.ctt.scm_numerical.scm_view import SCMView
from os.path import join
from matresdev.db import SimDB
from matresdev.db.exdb import ExRun
from mathkit.mfn.mfn_line.mfn_line import MFnLineArray
from scipy.optimize import fsolve, fmin

if __name__ == '__main__':

    simdb = SimDB()

    # specify the path to the data file.
    SH4V1_6cm = join(simdb.exdata_dir,
                     'tensile_tests',
                     'dog_bone',
                     '2012-04-12_TT-12c-6cm-0-TU_SH4',
                     'TT-12c-6cm-0-TU-SH4-V1.DAT')

    SH4V2_6cm = join(simdb.exdata_dir,
                     'tensile_tests',
                     'dog_bone',
                     '2012-04-12_TT-12c-6cm-0-TU_SH4',
                     'TT-12c-6cm-0-TU-SH4-V2.DAT')

    SH4V3_6cm = join(simdb.exdata_dir,
                     'tensile_tests',
                     'dog_bone',
                     '2012-04-12_TT-12c-6cm-0-TU_SH4',
                     'TT-12c-6cm-0-TU-SH4-V3.DAT')
    SH3V1_6cm = join(simdb.exdata_dir,
                     'tensile_tests',
                     'dog_bone',
                     '2012-03-20_TT-12c-6cm-0-TU_SH3',
                     'TT-12c-6cm-TU-0-SH3-V1.DAT')

    SH3V2_6cm = join(simdb.exdata_dir,
                     'tensile_tests',
                     'dog_bone',
                     '2012-03-20_TT-12c-6cm-0-TU_SH3',
                     'TT-12c-6cm-TU-0-SH3-V2.DAT')

    SH3V3_6cm = join(simdb.exdata_dir,
                     'tensile_tests',
                     'dog_bone',
                     '2012-03-20_TT-12c-6cm-0-TU_SH3',
                     'TT-12c-6cm-TU-0-SH3-V3.DAT')

    SH3V1_4cm = join(simdb.exdata_dir,
                     'tensile_tests',
                     'dog_bone',
                     '2012-03-20_TT-12c-4cm-0-TU_SH3',
                     'TT-12c-4cm-TU-0-SH3-V1.DAT')
    SH3V2_4cm = join(simdb.exdata_dir,
                     'tensile_tests',
                     'dog_bone',
                     '2012-03-20_TT-12c-4cm-0-TU_SH3',
                     'TT-12c-4cm-TU-0-SH3-V2.DAT')
    SH3V3_4cm = join(simdb.exdata_dir,
                     'tensile_tests',
                     'dog_bone',
                     '2012-03-20_TT-12c-4cm-0-TU_SH3',
                     'TT-12c-4cm-TU-0-SH3-V3.DAT')

    # construct the experiment
    V14 = ExRun(data_file=SH4V1_6cm)
    V24 = ExRun(data_file=SH4V2_6cm)
    V34 = ExRun(data_file=SH4V3_6cm)
    V13 = ExRun(data_file=SH3V1_6cm)
    V23 = ExRun(data_file=SH3V2_6cm)
    V33 = ExRun(data_file=SH3V3_6cm)
    V1small = ExRun(data_file=SH3V1_4cm)
    V2small = ExRun(data_file=SH3V2_4cm)
    V3small = ExRun(data_file=SH3V3_4cm)

    class Opt():

        r = 0.00345
        l = 1.0
        theta = 0.0
        xi = .013
        phi = 1.
        E_m = 25e3
        E_f = 200e3
        length = 1500.
        nx = 1500
        smallest_resid = 300.
        t_shape = 1.
        mtrx_scale = 5.8
        mtrx_shape = 21.0

        def get_eps_extrap(self, t_scale, ef):
            tau = RV('weibull_min', shape=self.t_shape, scale=t_scale)
            random_field = RandomField(seed=True,
                                       lacor=5.,
                                        xgrid=np.linspace(0., self.length, 500),
                                        nsim=1,
                                        loc=.0,
                                        shape=self.mtrx_shape,
                                        scale=self.mtrx_scale * 1.1,
                                        distribution='Weibull'
                                       )
            rf = CBEMClampedFiberStressSP()
            rand = FunctionRandomization(q=rf,
                                         tvars=dict(tau=tau,
                                                    l=self.l,
                                                    E_f=ef,
                                                    theta=self.theta,
                                                    xi=self.xi,
                                                    phi=self.phi,
                                                    E_m=self.E_m,
                                                    r=self.r,
                                                    V_f=0.0166
                                                         ),
                                            n_int=500
                                            )

            scm = SCM(length=self.length,
                      nx=self.nx,
                      random_field=random_field,
                      cb_randomization=rand,
                      cb_type='mean',
                      load_sigma_c_min=.1,
                      load_sigma_c_max=20.,
                      load_n_sigma_c=100,
                      n_w=50,
                      n_x=61,
                      n_BC=2
                      )

            scm_view = SCMView(model=scm)
            scm_view.model.evaluate()
            return scm_view.eps_sigma

        def get_eps(self, t_scale, ef):
            tau = RV('weibull_min', shape=self.t_shape, scale=t_scale)
            random_field = RandomField(seed=True,
                                       lacor=5.,
                                        xgrid=np.linspace(0., self.length, 500),
                                        nsim=1,
                                        loc=.0,
                                        shape=self.mtrx_shape,
                                        scale=self.mtrx_scale,
                                        non_negative_check=True,
                                        distribution='Weibull'
                                       )
            rf = CBEMClampedFiberStressSP()
            rand = FunctionRandomization(q=rf,
                                         tvars=dict(tau=tau,
                                                    l=self.l,
                                                    E_f=ef,
                                                    theta=self.theta,
                                                    xi=self.xi,
                                                    phi=self.phi,
                                                    E_m=self.E_m,
                                                    r=self.r,
                                                    V_f=0.0111
                                                         ),
                                            n_int=500
                                            )

            scm = SCM(length=self.length,
                      nx=self.nx,
                      random_field=random_field,
                      cb_randomization=rand,
                      cb_type='mean',
                      load_sigma_c_min=.1,
                      load_sigma_c_max=12.,
                      load_n_sigma_c=100,
                      n_w=50,
                      n_x=61,
                      n_BC=2
                      )

            scm_view = SCMView(model=scm)
            scm_view.model.evaluate()
            return scm_view.eps_sigma

        def experiment(self, e):
            eps = V13.ex_type.eps_asc[108:]
            sig = V13.ex_type.sig_c_asc[108:]
            line = MFnLineArray(xdata=eps, ydata=sig)
            vect_line = np.vectorize(line.get_value)
            return vect_line(e)

        def residuum(self, params):
            t_scale = params[0]
            ef = params[1]
            try:
                eps_opt, sigma = self.get_eps(t_scale, ef)
                squares = (self.experiment(eps_opt) - sigma) ** 2
                weight = np.ones(len(squares))
                weight[50:70] = 5.
                resid = np.sum((self.experiment(eps_opt) - sigma) ** 2 * weight)
            except:
                print '!!!!FAILED!!!!'
                resid = self.smallest_resid * 20.0
            print 't_scale', t_scale
            print 'Ef', ef
            print 'RESIDUUM', resid
            if resid < self.smallest_resid:
                self.smallest_resid = resid
                plt.figure()
                plt.plot(V13.ex_type.eps_asc, V13.ex_type.sig_c_asc, color='black')
                plt.plot(V23.ex_type.eps_asc, V23.ex_type.sig_c_asc, color='black')
                plt.plot(V33.ex_type.eps_asc, V33.ex_type.sig_c_asc, color='black')
                plt.plot(eps_opt, sigma, color='red',
                         label='%.2f, %.2f' %(t_scale, ef))
                plt.plot(V2small.ex_type.eps_asc, V2small.ex_type.sig_c_asc, color='black')
                plt.plot(V3small.ex_type.eps_asc, V3small.ex_type.sig_c_asc, color='black')
                try:
                    eps_extrap, sigma_extrap = self.get_eps_extrap(t_scale,
                                                                   ef)
                    plt.plot(eps_extrap, sigma_extrap, color='red')
                except:
                    pass
                plt.xlim(0., 0.008)
                plt.ion()
                plt.legend(loc='best')
            plt.show()
            return resid

        def get_parameters(self):
            opt_params = fmin(self.residuum, np.array([0.087, 201e3]))
            print opt_params

    def plot():

        # access the response values.
        for V in [V13, V23, V33]:
            eps = V.ex_type.eps_asc * 100.
            sig_c = V.ex_type.sig_c_asc
            plt.plot(eps, sig_c, color='black')
#        for V in [V2small, V3small]:
#            eps = V.ex_type.eps_asc * 100.
#            sig_c = V.ex_type.sig_c_asc
#            plt.plot(eps, sig_c, color='black')
        plt.xlim(0.0, 0.8)
        plt.ylim(0.0,27.)

    # filaments
        r = 0.00345
        tau = RV('weibull_min', shape=1.0, scale=.086)
        Ef = 195e3
        Em = 25e3
        l = 1.0#RV('uniform', scale=6., loc=.5)
        theta = 0.0
        xi = 20.013#RV('weibull_min', scale=0.02, shape=5)
        phi = 1.

        length = 1500.
        nx = 1500
        random_field1 = RandomField(seed=True,
                                   lacor=10.,
                                    xgrid=np.linspace(0., length, 600),
                                    nsim=1,
                                    loc=.0,
                                    shape=21.,
                                    scale=5.8 * 1.05,
                                    non_negative_check=True,
                                    distribution='Weibull'
                                   )

        rf = CBEMClampedFiberStressSP()
        rand1 = FunctionRandomization(q=rf,
                                     tvars=dict(tau=tau,
                                                l=l,
                                                E_f=Ef,
                                                theta=theta,
                                                xi=xi,
                                                phi=phi,
                                                E_m=Em,
                                                r=r,
                                                V_f=0.0166
                                                     ),
                                        n_int=20
                                        )

        scm1 = SCM(length=length,
                  nx=nx,
                  random_field=random_field1,
                  cb_randomization=rand1,
                  cb_type='mean',
                  load_sigma_c_min=.1,
                  load_sigma_c_max=25.,
                  load_n_sigma_c=200,
                  n_w=60,
                  n_x=81,
                  n_BC=2
                  )

        scm_view1 = SCMView(model=scm1)
        scm_view1.model.evaluate()

        random_field2 = RandomField(seed=True,
                                   lacor=10.,
                                    xgrid=np.linspace(0., length, 600),
                                    nsim=1,
                                    loc=.0,
                                    shape=21.,
                                    scale=5.8,
                                    non_negative_check=True,
                                    distribution='Weibull'
                                   )

        rand2 = FunctionRandomization(q=rf,
                                     tvars=dict(tau=tau,
                                                l=l,
                                                E_f=Ef,
                                                theta=theta,
                                                xi=xi,
                                                phi=phi,
                                                E_m=Em,
                                                r=r,
                                                V_f=0.0111
                                                     ),
                                        n_int=20
                                        )

        scm2 = SCM(length=length,
                  nx=nx,
                  random_field=random_field2,
                  cb_randomization=rand2,
                  cb_type='mean',
                  load_sigma_c_min=.1,
                  load_sigma_c_max=15.,
                  load_n_sigma_c=200,
                  n_w=60,
                  n_x=81,
                  n_BC=2
                  )

        scm_view2 = SCMView(model=scm2)
        scm_view2.model.evaluate()

        eps1, sigma1 = scm_view1.eps_sigma
        plt.plot(eps1 * 100., sigma1, lw=2, ls='dashed', color='red')
                 #label='no of cracks = ' + str(len(scm_view1.crack_widths(12.))))
        plt.legend(loc='best')
#        plt.xlabel('composite strain [-]')
#        plt.ylabel('composite stress [MPa]')
        eps2, sigma2 = scm_view2.eps_sigma
        plt.plot(eps2 * 100., sigma2, lw=2, color='red')
                 #label='no of cracks = ' + str(len(scm_view2.crack_widths(12.))))
        plt.legend(loc='best')
        plt.xlabel('composite strain [%]')
        plt.ylabel('composite stress [MPa]')
#        plt.legend(loc='best')
        #plt.figure()
        #plt.hist(scm_view.crack_widths(15.), bins=20, label='load = 15 MPa')
        #plt.hist(scm_view1.crack_widths(22.), bins=20, label='load = 22 MPa')
        #plt.legend(loc='best')
        plt.show()

    plot()
    #opt = Opt()
    #opt.get_parameters()
