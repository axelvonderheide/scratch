from stats.spirrid.spirrid import SPIRRID
from matplotlib import pyplot as plt
from stats.spirrid.rf_filament import Filament
from math import pi, factorial
from itertools import combinations, chain

def choose( n, k ):
    '''
        Combination coefficient
    '''
    return factorial( n ) / ( factorial( k ) * factorial( n - k ) )

def powerset( iterable ):
    '''
        Return object of all combination of iterable. 
        powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)
    '''
    s = list( iterable )
    return chain.from_iterable( combinations( s, r ) for r in range( len( s ) + 1 ) )

def run():
    # Quantities for the response function
    # and randomization
    # 
    #E_mod = 70 * 1e+9 # Pa
    #sig_u = 1.25 * 1e+9 # Pa
    D = 26 * 1.0e-6 # m
    A = ( D / 2.0 ) ** 2 * pi
    #xi_u = sig_u / E_mod

    dict_var = {'E_mod':70.e9,
                'xi':0.02,
                'A':A,
                'theta':0,
                'lambd':0}

    # random variable dictionary
    n_int = 300
    dict_rv = {
               'xi':{'variable':'xi', 'distribution':'weibull_min', 'scale':0.02, 'shape':10., 'n_int':n_int},
               'E_mod':{'variable':'E_mod', 'distribution':'uniform', 'loc':70e+9, 'scale':15e+9, 'n_int':n_int},
               'theta':{'variable':'theta', 'distribution':'uniform', 'loc':0.0, 'scale':0.01, 'n_int':n_int},
               'lambd':{'variable':'lambd', 'distribution':'uniform', 'loc':0.0, 'scale':.2, 'n_int':n_int},
               'A':{'variable':'A', 'distribution':'uniform', 'loc':A * 0.3, 'scale':0.7 * A, 'n_int':n_int},
               }

    # list of all combinations of response function parameters
    rv_comb_list = list( powerset( dict_rv ) )

    # define a tables with the run configurations to start in a batch
    run_list = [
#                ( 
#                 {'cached_dG'         : False,
#                  'compiled_QdG_loop'  : True,
#                  'compiled_eps_loop' : True },
#                  'bx-',
#                  '$\mathrm{C}_{e,\\theta} ( q(e,\\theta) \cdot g[\\theta_1]  g[\\theta_2] \dots g[\\theta_n] ) $ - %4.2f sec'
#                 ),
#                ( 
#                 {'cached_dG'         : False,
#                  'compiled_QdG_loop'  : True,
#                  'compiled_eps_loop' : False },
#                 'r-2',
#                 '$\mathrm{P}_{e} ( \mathrm{C}_{\\theta} ( q(e,\\theta) \cdot g[\\theta_1]  g[\\theta_2] \dots g[\\theta_n] ) ) $ - %4.2f sec',
#                 ),
                 ( 
                 {'cached_dG'         : True,
                  'compiled_QdG_loop'  : False,
                  'compiled_eps_loop' : False },
                 'y--',
                 '$\mathrm{P}_{e} ( \mathrm{N}_{\\theta} ( q(e,\\theta) \cdot G[\\theta] ) ) $ - %4.2f sec'
                 ),
#                ( 
#                 {'cached_dG'         : True,
#                  'compiled_QdG_loop'  : True,
#                  'compiled_eps_loop' : True },
#                  'go-',
#                  '$\mathrm{C}_{e,\\theta} ( q(e,\\theta) \cdot G[\\theta] ) $ - %4.2f sec',
#                  ),
                ]

    outfile = open( 'filament_comb_stat.dat', 'w' )
    sp1_outfile = open( 'filament_comb_1.dat', 'w' )
    sp2_outfile = open( 'filament_comb_2.dat', 'w' )

    plot = False

    for id, rv_comb in enumerate( rv_comb_list[15:25] ):#[0,-1]
        for i in range( 0, 1 ):#2
            num_rv = len( rv_comb ) # number of randomized variables
            if i == 0:
                sp1_outfile.write( '%i \t %i' % ( id, num_rv ) )
            else:
                sp2_outfile.write( '%i \t %i' % ( id, num_rv ) )
            # construct a default response function
            rf = Filament()
            rf.set( **dict_var )

            # construct the integrator and provide it with the response function.  
            s = SPIRRID( rf = rf,
                        min_eps = 0.00, max_eps = 0.04, n_eps = 20 ) # max_eps = 0.05

            for rv in rv_comb:
                s.add_rv( **dict_rv[rv] )

            if plot == True:
                legend = []

            print 'Combination', rv_comb, ' -- ', i
            outfile.write( "Combination %s -- %i\n" % ( rv_comb, i ) )
            outfile.write( "%s \t %s\n" % ( "Run", "Time" ) )

            for idx, run in enumerate( run_list ):
                run_options, plot_options, legend_string = run
                print 'run', idx,
                s.set( **run_options )
                if plot == True:
                    plt.figure( id )
                    s.mean_curve.plot( plt, plot_options )
                    legend.append( legend_string % s.exec_time )
                print 'execution time', s.exec_time
                outfile.write( "%s \t %.6f\n" % ( idx, s.exec_time ) )

                if i == 0:
                    sp1_outfile.write( "\t %.6f" % ( s.exec_time ) )
                else:
                    sp2_outfile.write( "\t %.6f" % ( s.exec_time ) )

            outfile.write( "=================\n" )
            outfile.flush()

            if i == 0:
                sp1_outfile.write( '\t %s\n' % str( rv_comb ) )
                sp1_outfile.flush()
            else:
                sp2_outfile.write( '\t %s\n' % str( rv_comb ) )
                sp2_outfile.flush()
            if plot == True:
                plt.xlabel( 'strain [-]' )
                plt.ylabel( 'stress' )
                plt.legend( legend )
                plt.title( s.rf.title )

        del s


    outfile.close()
    sp1_outfile.close()
    sp2_outfile.close()
    if plot == True:
        plt.show()

if __name__ == '__main__':
    run()










