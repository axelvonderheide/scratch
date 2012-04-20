
from stats.spirrid import SPIRRID
from matplotlib import pyplot as plt
from rf_simple import Filament
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
    dict_var = {'x' : 10,
                'y' : 20,
                'z' : 30,
                'a' : 40}
    
    n_var = len( dict_var )
    n_com = 0
    for k in range( 0, n_var ):
        n_com += choose( n_var, k )
        
    outfile = open( 'simple_comb_stat.dat', 'w' )
    
    for comb in range( 1, n_com + 1 ):
        for i in range( 0, 2 ):
            # construct a default response function for a single filament
            rf = Filament()
            rf.set( **dict_var ) # (x=10, y=20, ...)
        
            # construct the integrator and provide it with the response function.
        
            s = SPIRRID( rf=rf,
                         min_eps=0.00, max_eps=0.05, n_eps=1 )
        
            # construct the random variables
        
            n_int = 35
            
            dict_rv = {
                       'x':{'variable':'x', 'distribution':'norm', 'loc':10., 'scale':2., 'n_int':n_int},
                       'y':{'variable':'y', 'distribution':'norm', 'loc':10., 'scale':2., 'n_int':n_int},
                       'z':{'variable':'z', 'distribution':'norm', 'loc':10., 'scale':2., 'n_int':n_int},
                       'a':{'variable':'a', 'distribution':'norm', 'loc':10., 'scale':2., 'n_int':n_int},
                       }
            list_rv_comb = list( powerset( dict_rv ) )
            for rv in list_rv_comb[comb]:
                s.add_rv( **dict_rv[rv] )
                    
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
        
            print 'Combination', list_rv_comb[comb], ' -- ', i
            outfile.write( "Combination %s -- %i\n" % ( list_rv_comb[comb], i ) )
            outfile.write( "%s \t %s\n" % ( "Run", "Time" ) )
            
            for idx, run in enumerate( run_list ):
                run_options, plot_options, legend_string = run
                print 'run', idx,
                s.set( **run_options )
                s.mean_curve.plot( plt, plot_options )
                print 'execution time', s.exec_time, s.mu_q_peak
                outfile.write( "%s \t %s\n" % ( idx, s.exec_time ) )
                legend.append( legend_string % s.exec_time )

            
            outfile.write( "=================\n" )

            plt.xlabel( 'strain [-]' )
            plt.ylabel( 'stress' )
            plt.legend( legend )
        
            plt.title( s.rf.title )
            del s 
    plt.show()
    outfile.close()

if __name__ == '__main__':
    run()

