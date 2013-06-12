import numpy as np
from matplotlib import pyplot as plt
from stats.spirrid.rv import RV
import math
import numpy as np

###################
def same_length( l ):
    return np.arccos( l[0] / l )
    
n_int = 200

phi = RV( 'sin2x', loc = 0., scale = 1. )
a_l = RV( 'uniform', loc = 0., scale = 1. )
discr_ppf = np.linspace( .0001 / n_int, 1. - .0001 / n_int, n_int )
phi = phi.ppf( discr_ppf )
a = a_l.ppf( discr_ppf )
#for ae in a:
#    x = np.ones( len( a ) ) * ae
#    print x, phi
#    plt.plot( x , phi, 'ro' )
    
isol = np.linspace( a[5], a[-5], 5 )
'''
for i, il in enumerate( isol ):
    
    x_l = np.linspace( il, np.max( a ), 200 )
    y_l = same_length( x_l )
    if i != 1:
        plt.plot( x_l, y_l , 'k' )
    else:
        x_f = x_l
        y_f = y_l
        plt.plot( x_l, y_l , 'r' )
        '''
        
#plt.plot( a, phi )
plt.show()
###########################################
'''Graph mit Erklaerungen'''
##########################################
idx = 0
for i, il in enumerate( isol ):
    
    x_l = np.linspace( il, np.max( a ), 200 )
    y_l = same_length( x_l )
    if i != 1:
        pass
        #plt.plot( x_l, y_l , color = '0.01' , linewidth = '0.5' )
    else:
        x_f = x_l
        y_f = y_l
        idx = ( np.abs( y_l - 2.5 * math.pi / 8 ) ).argmin()
        plt.axvline( x_f[idx], ymin = 0, ymax = ( 0.99 * math.pi / 5. ), linewidth = '0.5' , color = 'k' , linestyle = '-' )
        plt.plot( x_l, y_l , 'k', linewidth = '2' )
        plt.axhline( y = 2.5 * math.pi / 8., xmin = 0, xmax = 1, color = 'k', linewidth = '2' , linestyle = '--' )
        plt.yticks( ( 0, math.pi * 2.5 / 8., math.pi / 2. ),
               ( "0", r"$\phi(\chi_f)$", r"$\frac{\pi}{2}$" ), fontsize = '20' )
        plt.xticks( ( x_f[0], x_f[idx], 1 ), ( "$a$", "$a_s$", r"$\i_f$" ), fontsize = '20' )
        plt.xlim( 0, 1 )
        #plt.grid( True )
        plt.ylim( 0, math.pi / 2. )
        plt.ylabel( r"$\phi$", fontsize = '24' )
        plt.xlabel( r"$a_e$", fontsize = '24' )
        x_fill_list = list( x_f[:idx] )
        x_fill_list.append( x_f[idx] + 1e-15 )
        y_fill_list = list( y_f[:idx] )
        y_fill_list.append( 0 )
        plt.fill( x_fill_list, y_fill_list , hatch = '/' , fill = False, linewidth = 0.0001 )
        x_fillB_list = list( x_f[idx:] )
        x_fillB_list.append( x_f[-1] + 1e-15 )
        x_fillB_list.insert( 0, x_f[idx] )
        y_fillB_list = list( np.ones( len( y_f[idx:] ) ) * y_f[idx] )
        y_fillB_list.append( 0 )
        y_fillB_list.insert( 0, 0 )
        plt.fill( x_fillB_list, y_fillB_list , hatch = ' \ ' , fill = False, linewidth = 0.0001 )
    
    bbox_props = dict( boxstyle = "round", fc = "w", ec = "0.5", alpha = 0.9 )
    plt.text( 0.39, 0.5, "A", ha = "center", va = "center", size = 30,
        bbox = bbox_props )
    plt.text( 0.73, 0.5, "B", ha = "center", va = "center", size = 30,
        bbox = bbox_props )
    
        

plt.show()
