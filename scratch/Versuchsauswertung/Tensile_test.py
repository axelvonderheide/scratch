'''
Created on 03.07.2013

@author: acki
'''
import numpy as np
from scipy.interpolate import interp1d
from quaducom.micro.resp_func.CB_rigid_mtrx import CBResidual
import numpy as np
from spirrid.spirrid import SPIRRID
from spirrid.rv import RV
from matplotlib import pyplot as plt
from scipy.optimize import minimize
def TensileTest( Dataplot, meanplot ):
    l_of_all_data = ['DATA/PO01_RYP.ASC', 'DATA/PO02_RYP.ASC', 'DATA/PO03_RYP.ASC', \
              'DATA/PO04_RYP.ASC', 'DATA/PO05_RYP.ASC', 'DATA/PO06_RYP.ASC', 'DATA/PO07_RYP.ASC', \
              'DATA/PO08_RYP.ASC', 'DATA/PO09_RYP.ASC']
    
    w_x = np.linspace( 0, 18, 300 )
    # use to enter in Data_to_plot if you need all
    all_carbon = [1, 2, 3, 4, 8]
    all_glass = [5, 6, 7, 9]
    ####################################
    '''Enter Pullouts u want to plot'''
    ####################################
    Data_to_plot = np.array( [1, 2, 3, 4] )
    
    ####################################
    
    # Data into arrays
    plotlist = []
    for i in Data_to_plot:
        plotlist.append( l_of_all_data[i - 1] )
    mean = []
    F_w_list = []
    for i, dataname in enumerate( plotlist ):
        f = open( dataname, 'r' )
        time = []
        F = []
        w1 = []
        w2 = []
        for line in f:
                            #1###
                            index1 = line.index( ';' ) 
                            time.append( float( line[:index1] ) )
                            line = line[index1 + 1:]
                            #2###
                            index2 = line.index( ';' )
                            F.append( np.float( line[:index2] ) )
                            line = line[index2 + 1:]
                            #3###
                            index3 = line.index( ';' )
                            w1.append( np.float( line[:index3] ) )
                            line = line[index3 + 1:]
                            #4###
                            w2.append( np.float( line ) )
                            
        F = np.array( F )
        w1 = np.array( w1 )
        w2 = np.array( w2 )
        F = F - F[0] + 0.05
        w2 = -w2[0] - w2
        w1 = -w1[0] - w1
        w_mid = ( w1 + w2 ) / 2
        # plt.plot( w_mid, F )
        string = '%d'.format( i )
        # plt.plot( w2, F )
        f.close()
        F_w_item = interp1d( w_mid, F, bounds_error = False, fill_value = 0 )
        F_w_list.append( F_w_item )
        if Dataplot == True:
            plt.plot( w_x, F_w_item( w_x ), label = 'Pullout: {}'.format( Data_to_plot[i] ), linewidth = '2' )
        if len( mean ) == 0:
            mean = F_w_item( w_x )
        else:
            mean += F_w_item( w_x )
    if meanplot == True: 
        mean = mean / len( F_w_list )  
        plt.plot( w_x, mean, label = 'mean', linewidth = '3', color = 'k' )
    
    
    ####################################MODEL##########################################
    def Model( itertuple ):
        tau_shape, tau_scale = itertuple
        def CB_composite_stress( w, tau, E_f, V_f, r, m, sV0, Pf, n_int ):
            cb = CBResidual()
            spirrid = SPIRRID( q = cb,
                        sampling_type = 'PGrid',
                        eps_vars = dict( w = w ),
                        theta_vars = dict( tau = tau, E_f = E_f, V_f = V_f, r = r,
                                   m = m, sV0 = sV0, Pf = Pf ),
                        n_int = n_int )
            if isinstance( r, RV ):
                r_arr = np.linspace( r.ppf( 0.001 ), r.ppf( 0.999 ), 300 )
                Er = np.trapz( r_arr ** 2 * r.pdf( r_arr ), r_arr )
            else:
                Er = r ** 2
            sigma_c = spirrid.mu_q_arr / Er    
            # plt.ylim(0, 35)
            return w, sigma_c * 0.445 / 1e3 * 1300 
        w = np.linspace( 0, 18, 300 )
        tau = RV( 'weibull_min', shape = tau_shape, scale = tau_scale )
        E_f = 170e3
        V_f = .001
        r = 0.03  # RV('uniform', loc=0.001, scale=0.004)
        m = 7.
        # sV0=XXX corresponds to sL0=0.02 at L0=100 and r=0.002
        sV0 = 3.0e-3
        Pf = RV( 'uniform', loc = 0., scale = 1.0 )
        n_int = 30
    
        return CB_composite_stress( w, tau, E_f, V_f, r, m, sV0, Pf, n_int )
    count = 0
    weighted_array = np.ones( len( w_x ) ) * ( w_x < 5 ) + ( w_x >= 5 ) * 0.0001
    def residuum( itertuple ): 
        print itertuple
        # print mean - Model( itertuple )
        return np.sum( ( mean - Model( itertuple ) ) ** 2 * weighted_array )
    
    # iterF = minimize( residuum, x0 = [3., 0.03], method = 'nelder-mead' )
    # t_shape, t_scale = iterF
    

    
    w, sigma_c = Model( [  3.59452024e+01 , 4.89165913e-03] )
    plt.plot( w, sigma_c )
    plt.legend()
    plt.show()
    
    
#################Plots
    
Dataplot = False
meanplot = True
    
    
    
    

TensileTest( Dataplot, meanplot )
