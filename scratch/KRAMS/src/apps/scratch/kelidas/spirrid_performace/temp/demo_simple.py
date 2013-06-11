
from stats.spirrid import SPIRRID
from matplotlib import pyplot as plt
from rf_simple import Filament
from math import pi

def run():
    # Quantities for the response function
    # and randomization
    # 
    xx = 10 # Pa
    yy = 20 # Pa
    zz = 30
    aa = 40


    # construct a default response function for a single filament

    rf = Filament( x=xx, y=yy, z=zz, a=aa )

    # construct the integrator and provide it with the response function.

    s = SPIRRID( rf=rf,
                 min_eps=0.00, max_eps=0.05, n_eps=1 )

    # construct the random variables

    n_int = 80

    s.add_rv( 'x', distribution='norm', scale=10., shape=2., n_int=n_int )
    s.add_rv( 'y', distribution='norm', loc=10., scale=2., n_int=n_int )
    s.add_rv( 'z', distribution='norm', loc=10., scale=2., n_int=n_int )
    s.add_rv( 'a', distribution='norm', loc=10., scale=2., n_int=n_int )

    
    # define a tables with the run configurations to start in a batch

    run_list = [
                ( 
                 {'cached_qg'         : False,
                  'compiled_qg_loop'  : True,
                  'compiled_eps_loop' : True },
                  'bx-',
                  '$\mathrm{C}_{e,\\theta} ( q(e,\\theta) \cdot g[\\theta_1]  g[\\theta_2] \dots g[\\theta_n] ) $ - %4.2f sec'
                 ),
                ( 
                 {'cached_qg'         : False,
                  'compiled_qg_loop'  : True,
                  'compiled_eps_loop' : False },
                 'r-2',
                 '$\mathrm{P}_{e} ( \mathrm{C}_{\\theta} ( q(e,\\theta) \cdot g[\\theta_1]  g[\\theta_2] \dots g[\\theta_n] ) ) $ - %4.2f sec',
                 ),
                ( 
                 {'cached_qg'         : True,
                  'compiled_qg_loop'  : True,
                  'compiled_eps_loop' : True },
                  'go-',
                  '$\mathrm{C}_{e,\\theta} ( q(e,\\theta) \cdot G[\\theta] ) $ - %4.2f sec',
                  ),
                ( 
                 {'cached_qg'         : True,
                  'compiled_qg_loop'  : False,
                  'compiled_eps_loop' : False },
                 'y--',
                 '$\mathrm{P}_{e} ( \mathrm{N}_{\\theta} ( q(e,\\theta) \cdot G[\\theta] ) ) $ - %4.2f sec'
                 ),
                ]

    legend = []

    for idx, run in enumerate( run_list ):
        run_options, plot_options, legend_string = run
        print 'run', idx,
        s.set( **run_options )
        s.mean_curve.plot( plt, plot_options )
        print 'execution time', s.exec_time
        legend.append( legend_string % s.exec_time )

    plt.xlabel( 'strain [-]' )
    plt.ylabel( 'stress' )
    plt.legend( legend )

    plt.title( s.rf.title )
    plt.show()

if __name__ == '__main__':
    run()
