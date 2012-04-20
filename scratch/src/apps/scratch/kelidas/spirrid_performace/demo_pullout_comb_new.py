from stats.spirrid import \
    SPIRRID, RV
from matplotlib import pyplot as plt
from quaducom.old.pullout.constant_friction_finite_fiber import ConstantFrictionFiniteFiber
from math import pi, factorial
from itertools import combinations, chain
import numpy as np

def choose(n, k):
    '''
        Combination coefficient
    '''
    return factorial(n) / (factorial(k) * factorial(n - k))

def powerset(iterable):
    '''
        Return object of all combination of iterable. 
        powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)
    '''
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s) + 1))

def run():
    # Quantities for the response function
    # and randomization
    # 
    #E_mod = 70 * 1e+9 # Pa
    #sig_u = 1.25 * 1e+9 # Pa
    #D = 26 * 1.0e-6 # m
    #A = ( D / 2.0 ) ** 2 * pi
    #xi_u = sig_u / E_mod

    # construct a default response function for a single filament
    dict_var = {'fu':1200.0e6,
                'qf':1500.,
                'L':0.02,
                'A':5.30929158457e-10,
                'E_mod':70.0e9,
                'z':0.0,
                'phi':0.0,
                'f':0.0}

    e_arr = np.linspace(0, 0.012, 40)

    # random variable dictionary
    tvars = dict(fu = RV('weibull_min', 1200.0e6, 200.),
                   qf = RV('uniform', 1500., 100.),
                   L = RV('uniform', 0.02, 0.02 / 2.),
                  A = RV('uniform', 5.30929158457e-10, .03 * 5.30929158457e-10),
                  E_mod = RV('uniform', 70.e9, 250.e9),
                  z = RV('uniform', 0.0, 0.03),
                   phi = RV('sin_distr', 0.0, 1.0),
                  f = RV('uniform', 0.0, 0.03))

    # list of all combinations of response function parameters
    rv_comb_list = list(powerset(dict_var))

    # define a tables with the run configurations to start in a batch
    run_list = [
                ('c',
                 {'cached_dG'         : False,
                  'compiled_eps_loop' : True },
                  'bx-',
                  '$\mathrm{C}_{e,\\theta} ( q(e,\\theta) \cdot g[\\theta_1]  g[\\theta_2] \dots g[\\theta_n] ) $ - %4.2f sec'
                 ),
                ('c',
                 {'cached_dG'         : False,
                  'compiled_eps_loop' : False },
                 'r-2',
                 '$\mathrm{P}_{e} ( \mathrm{C}_{\\theta} ( q(e,\\theta) \cdot g[\\theta_1]  g[\\theta_2] \dots g[\\theta_n] ) ) $ - %4.2f sec',
                 ),
                ('numpy',
                 { },
                  'go-',
                  '$\mathrm{C}_{e,\\theta} ( q(e,\\theta) \cdot G[\\theta] ) $ - %4.2f sec',
                  ),
                 ('c',
                 {'cached_dG'         : True,
                  'compiled_eps_loop' : False },
                 'y--',
                 '$\mathrm{P}_{e} ( \mathrm{N}_{\\theta} ( q(e,\\theta) \cdot G[\\theta] ) ) $ - %4.2f sec'
                 ),
                ]

    outfile = open('pullout_comb_stat_new.dat', 'w')
    sp1_outfile = open('pullout_comb_1_new.dat', 'w')
    sp2_outfile = open('pullout_comb_2_new.dat', 'w')

    plot = True

    for id, rv_comb in enumerate(rv_comb_list[163:219]):#[0,-1]
        for i in range(0, 2):
            num_rv = len(rv_comb) # number of randomized variables
            if i == 0:
                sp1_outfile.write('%i \t %i' % (id, num_rv))
            else:
                sp2_outfile.write('%i \t %i' % (id, num_rv))
            # construct a default response function

            # construct the integrator and provide it with the response function.  

            tvars_up = dict_var
            for rv in rv_comb:
                tvars_up[rv] = tvars[rv]

            if plot == True:
                legend = []

            s = SPIRRID(q = ConstantFrictionFiniteFiber(),
                e_arr = e_arr,
                n_int = 20,
                tvars = tvars_up,
                )
            s.sampling_type = 'TGrid'
            s.recalc = True
            s.sampling.recalc = True
            print 'Combination', rv_comb, ' -- ', i
            outfile.write("Combination %s -- %i\n" % (rv_comb, i))
            outfile.write("%s \t %s\n" % ("Run", "Time"))

            for idx, run in enumerate(run_list):
                code, run_options, plot_options, legend_string = run
                s.codegen_type = code
                s.codegen.set(**run_options)
                print 'run', idx,
                if plot == True:
                    plt.figure(id)
                    plt.plot(s.evar_list[0], s.mu_q_arr, plot_options)
                    legend.append(legend_string % s.exec_time)
                print 'execution time', s.exec_time
                outfile.write("%s \t %.6f\n" % (idx, s.exec_time))

                if i == 0:
                    sp1_outfile.write("\t %.6f" % (s.exec_time))
                else:
                    sp2_outfile.write("\t %.6f" % (s.exec_time))
            #plt.show()

            outfile.write("=================\n")
            outfile.flush()

            if i == 0:
                sp1_outfile.write('\t %s\n' % str(rv_comb))
                sp1_outfile.flush()
            else:
                sp2_outfile.write('\t %s\n' % str(rv_comb))
                sp2_outfile.flush()
            if plot == True:
                plt.xlabel('strain [-]')
                plt.ylabel('stress')
                plt.legend(legend)
                #plt.title(s.rf.title)


        del s


    outfile.close()
    sp1_outfile.close()
    sp2_outfile.close()
    if plot == True:
        plt.show()

if __name__ == '__main__':
    run()
