'''
Created on 4.1.2013

@author: Acki
'''
from rf_indep_CB_rand_xi import CBRandXi
from stats.spirrid.spirrid import SPIRRID
from stats.spirrid.rv import RV
from matplotlib import pyplot as plt
import numpy as np
from scipy.optimize import fmin
from etsproxy.mayavi import mlab
from stats.spirrid import make_ogrid as orthogonalize
from scipy.optimize import minimize
from math import pi
from scipy.special import gammainc, gamma
import pickle
def H( x ):
    return x >= 0.0

cb = CBRandXi()
 
def rand_tau_r( E_f, V_f, mu_tau, mu_r, m_list, CoV_tau_range, CoV_r_range, dp_tau, dp_r, sV0 ):
    #rand tau and rand r
    m = m_list[0]
    s0 = ( ( mu_tau * ( m + 1 ) * sV0 ** m ) / ( E_f * pi * mu_r ** 3 ) ) ** ( 1. / ( m + 1 ) )
    mu_xi = s0 * gamma( 1. + 1. / ( 1. + m ) )
    for i_m, mi in enumerate( m_list ):
        s0i = mu_xi / gamma( 1. + 1. / ( 1. + mi ) )
        sV0i = ( ( s0i ** ( mi + 1 ) * E_f * pi * mu_r ** 3. ) / ( mu_tau * ( mi + 1. ) ) ) ** ( 1. / mi )
        #
        #print 'Nr', i_m, 's0i', s0i, 'sV0i', sV0i, 'mu_xi', s0i * gamma(1. + 1. / (1. + mi))
        #
        #Pf = RV( 'uniform', loc = 0.0, scale = 1.0 )
        #w_arr = np.linspace(0, 1.2, 30)
        # loc scale generation for specific CoV
        #CoVtau
        CoV_tau_arr = np.linspace( CoV_tau_range[0], CoV_tau_range[1], dp_tau )
        loc_tau_arr = mu_tau - CoV_tau_arr * mu_tau * 3 ** 0.5
        scale_tau_arr = 2 * mu_tau - 2 * loc_tau_arr
        #CoVr
        CoV_r_arr = np.linspace( CoV_r_range[0], CoV_r_range[1], dp_r )
        loc_r_arr = mu_r - CoV_r_arr * mu_r * 3 ** 0.5
        scale_r_arr = 2 * mu_r - 2 * loc_r_arr
        #shaping for mayavi
        e_arr = orthogonalize( [CoV_tau_arr, CoV_r_arr] )
        x_axis = e_arr[0]
        y_axis = e_arr[1]
        #TAU gen Tuple of [loc,scale]
        stats_tau = []
        for s in range( dp_tau ):
            stats_tau.append( RV( 'uniform', loc = loc_tau_arr[s], scale = scale_tau_arr[s] ) )
        stats_tau[0] = mu_tau
        stats_r = []
        for s in range( dp_r ):
            stats_r.append( RV( 'uniform', loc = loc_r_arr[s], scale = scale_r_arr[s] ) )
        stats_r[0] = mu_r
        #r gen Tuple of [loc,scale]

        sigma_array = np.zeros( ( dp_tau, dp_r ) )
        w_array = np.zeros( ( dp_tau, dp_r ) )
        
        ##grid
        sigma_array_grid = np.zeros( ( dp_tau, dp_r ) )
        w_array_grid = np.zeros( ( dp_tau, dp_r ) )
        for i_tau, taui in enumerate( stats_tau ):
                for i_r, ri in enumerate( stats_r ):
                    total = SPIRRID( q = cb,
                        sampling_type = 'PGrid',
                        evars = dict(),
                        tvars = dict( tau = taui, E_f = E_f, V_f = V_f, r = ri,
                        m = mi, sV0 = sV0i ),
                        n_int = 70 )
                    
                    if isinstance( ri, RV ):
                        r_arr = np.linspace( ri.ppf( 0.001 ), ri.ppf( 0.999 ), 200 )
                        Er = np.trapz( r_arr ** 2 * ri.pdf( r_arr ), r_arr )
                    else:
                        Er = ri ** 2
                    #Optimization  taui, E_f, V_f, ri, mi, sV0i, Pf, Er
                    def spirrid_func( w_in, Er ):
                        w_arr = np.array( w_in )
                        #print 'w', w_arr
                        total.evars = dict( w = w_arr )
                        g = -1.*total.mu_q_arr / Er
                        #print 'res', g
                        total.evars = dict( w = w_arr )
                        return g
                    #w_test = np.linspace(0.001, 1, 1000)
                    #total.evars = dict(w=w_test)
                   # plt.plot(w_test, -1.*total.mu_q_arr / Er)
                    #plt.show()
                    print i_tau, 'von', dp_tau
                    #Optimierung
                    w_ult, sigma_ult, trash1, trash2, trash3 = fmin( spirrid_func, x0 = 0.1, args = [Er], ftol = 0.2, full_output = 1 )
                    sigma_array[i_tau, i_r] = np.float( -1.*sigma_ult )
                    w_array[i_tau, i_r] = np.float( w_ult )
                    
                    #grid
                    w_grid = np.linspace( 0.2, 0.28, 50 )
                    res_grid = -1.* spirrid_func( w_grid, Er )
                    index_max_grid = np.argmax( res_grid )
                    sigma_array_grid[i_tau, i_r] = res_grid[index_max_grid]
                    w_array_grid[i_tau, i_r] = w_grid[index_max_grid]
                    #print w_grid[index_max_grid]
                    
                   

        #print np.max( sigma_array )
        s_dataname = 'sigmaOPT802_with_m{}.npy'.format( mi )
        np.save( s_dataname, sigma_array )
        w_dataname = 'wOPT802_with_m{}.npy'.format( mi )
        np.save( w_dataname, w_array )
        
        #Grid Datasave
        s_gridname = 'sigmaGRID802_with_m{}.npy'.format( mi )
        np.save( s_gridname , sigma_array_grid )
        w_gridname = 'wGRID802_with_m{}.npy'.format( mi )
        np.save( w_gridname, w_array_grid )
        #mayaviplot
        mlab.surf( x_axis, y_axis, w_array )

    mlab.xlabel( "rand tau" )
    mlab.ylabel( "rand r" )
    mlab.zlabel( "sigma" )
    mlab.show()
    
###########################################################################################################


#Parameter


#PLOTS
#########
"Params"
#########
E_f, V_f, sV0 = 200e3, 0.01, 3.e-3

#########################
"random tau and random r"
#########################

#mu of rndm variables
mu_tau = 0.1
mu_r = 0.01

#plots of different m's with fix mu_xi
m_list = [3.]

#Ranges of CoV's
CoV_tau_range = [1e-7, 0.5]
CoV_r_range = [1e-7, 0.5]
dp_tau = 30
dp_r = 30

#Plot
rand_tau_r( E_f, V_f, mu_tau, mu_r, m_list, CoV_tau_range, CoV_r_range, dp_tau, dp_r, sV0 )
