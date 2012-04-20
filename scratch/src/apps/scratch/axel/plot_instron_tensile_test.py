'''
Created on 15.08.2011

@author: axel
'''
from scipy.optimize import curve_fit
from numpy import array, mean, var, sqrt, linspace, argmax, argmin, abs, sign
from math import e, pi as Pi
from quaducom.resp_func.cb_short_fiber import CBShortFiber
from matplotlib import pyplot as plt
import os


###############################################################################################################################################################
'''                            ''''''Hier Pfad angeben und die anzuzeigenden Versuchsreihen''''''                                                           '''
###############################################################################################################################################################

#Data Listcheck
path = "C:\\tensile_tests"  #Data path
dirList = os.listdir( path )


#specimen name
'''copy and paste'''
#name = 'pullout0', 'pullout15', 'pullout30'
#name = 'a', 'b', 'c', 'd', 'e', 'f'
name = ['d']


################################################################################################################################################################

def H( x ):
    return sign( sign( x ) + 1. )



#specimen number reaches to 20
n = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20']
stiffness_of_the_machine = 5500

colors = 'b', 'r', 'y', 'g', 'k', 'c', 'm'

weights0 = [6.9, 7.4, 7.4, 6.8, 6.9]
weights15 = [4, 4, 4, 4, 4, 4, 4]
weights30 = [6.4, 6.4, 7, 6.7]

cross_sections_a = [25.78 * 23.6, 27.5 * 27.1, 28.92 * 28.77, 27.2 * 26.6, 27.1 * 26.8, 29 * 26.25, 28.2 * 26.2, 28.1 * 28.94, 28.8 * 29.13, 28.2 * 26.0]
cross_sections_b = [28.34 * 23.96, 27.03 * 25.05, 23.76 * 23.42, 25.2 * 26.11, 25.13 * 26.73, 27.52 * 27.91, 26.84 * 26.57, 26.19 * 26.71, 25.5 * 26.3, 26.1 * 28.48]
cross_sections_c = [27.4 * 28.5, 26.53 * 25.56, 26.5 * 24, 29.6 * 25.3, 26.82 * 25.13, 26.73 * 26.71, 28.5 * 25.3, 30.4 * 26.5, 27.3 * 27.7, 25.66 * 25.25]
cross_sections_d = [25.12 * 27.17, 25.5 * 27.76, 26.77 * 25, 25.29 * 24.85, 27.45 * 29.2, 26.7 * 26.9, 25.94 * 26.73, 26.23 * 27.73, 30 * 30, 24 * 25.38]
cross_sections_e = [27 * 30, 25.5 * 27, 26.1 * 25.9, 24.5 * 26.2, 26.6 * 26.5, 29.5 * 27.5, 26 * 28.5]
cross_sections_f = [25.8 * 26.5, 26.1 * 26.9, 29.5 * 27.5, 26 * 28.5, 28 * 27.5, 26 * 28.5, 28 * 27.5, 31 * 25, 28.5 * 29 , 27 * 28]

###############################################

E = 200000.
A = 0.15 * 0.15

U = linspace( 0, 10, 100 )

#Datenfitting:
length = 100
width = 30
height = 30
L_f = 20
r = 0.075
Vf = 1

D_f = 0.15
E_f = 200000
le = 20
phi = 0

def fittingTAU( w, tau ):


        T = tau * Pi * D_f
        E = E_f
        A = D_f ** 2 / 4. * Pi
        #debonding stage
        q_deb = sqrt( E * A * T * w * H( w ) )
        #print q_deb
        # displacement at which debonding is finished
        w0 = le * T * ( le ) / E_f / A
        # pulling out stage - the fiber is pulled out from the
        # side with the shorter embedded length only
        q_pull = le * T * ( ( w0 - w ) / ( ( le + 1e-15 ) - w0 ) + 1 )
        q = q_deb * H( le * T - q_deb ) + q_pull * H( q_deb - le * T )

        # include inclination influence
        q = q * H( q ) #* e ** ( phi * f )


        return q


y_add = 0

###############################
color_run = 0
max_strength_list = []
fit_list_x = []
fit_list_y = []
tau_list = []
for specname in name:
    color = colors[color_run]
    color_run += 1
    for count in range( len( n ) )  :
        if specname + n[count] + '.raw' in dirList:
            f = open( '/tensile_tests/' + specname + n[count] + '.raw', 'r' )
            i = 0
            list_x = []
            list_y = []
            max_strength = 0
            for line in f:
                i += 1
                if i > 41:
                    #X-Data
                    tmp1 = line[line.index( ';' ) + 1:]
                    tmp1 = tmp1.replace( ',', '.' )
                    strain = tmp1[:tmp1.index( ';' )]
                    strain = strain.replace( ',', '.' )
                    list_x.append( float( strain ) )
                    #Y-Data
                    tmp2 = tmp1[tmp1.index( ';' ) + 1:]
                    tmp2 = tmp2[tmp2.index( ';' ) + 1:]
                    strength = float( tmp2[:tmp2.index( ';' )] )
                    list_y.append( strength )
            list_x = array( list_x )
            list_y = array ( list_y )
            #list_x = list_x - list_y / stiffness_of_the_machine
            if specname[0:7] == 'pullout':
                t_name = 'list_y=list_y-' + 'weights' + specname[7:] + '[' + str( count ) + ']'
                exec  t_name
            if len( specname ) <= 2:
                t_name = 'list_y=list_y/' + 'cross_sections_' + specname[0] + '[' + str( count ) + ']'
                exec t_name
                pass


            while list_y[0] < 0:
                index_zero = argmin( abs( list_y[:200] ) ) + 1
                list_y = list_y[index_zero:]
                list_x = list_x[index_zero:]
                if len( list_y ) == 0:
                    print specname, n[count]

            if list_x[0] < 0:
                list_x += abs( list_x[0] )
            else: list_x -= abs( list_x[0] )

            max_strength = max( list_y )
            max_strength_list.append( max_strength )

            fit_list_x = list_x[:argmax( list_y )]
            fit_list_y = list_y[:argmax( list_y )]
            plt.plot( list_x, list_y, color )



            if specname[0:7] == 'pullout':
                phi = float( specname[7:] )
                tau, vartau = curve_fit( fittingTAU, fit_list_x, fit_list_y )
                tau_list.append( tau )
    if specname[0:7] == 'pullout':
        #plt.plot( list_x, list_y, color  )
        print 'tau_mean', mean( array( tau_list ) ), 'tau_stdev', sqrt( var( array( tau_list ) ) ), 'tau', tau_list
    if specname[0] == 'a':
        plt.plot( list_x, list_y, color, label = '0% Vf' )
    if specname[0] == 'b':
        plt.plot( list_x, list_y, color, label = '2% Vf' )
    if specname[0] == 'c':
        plt.plot( list_x, list_y, color, label = '2.5% Vf' )
    if specname[0] == 'd':
        plt.plot( list_x, list_y, color, label = '3% Vf' )
    if specname[0] == 'e':
        plt.plot( list_x, list_y, color, label = '3.5% Vf' )
    if specname[0] == 'f':
        plt.plot( list_x, list_y, color, label = '4% Vf' )


    x = array( max_strength_list )
    print 'Versuchsreihe', specname, 'mean', mean( x ), 'Stdev', sqrt( var( x ) ) , 'max strengths', max_strength_list
    max_strength_list = []
    tau_list = []
    f.close()


plt.xlabel( 'Verschiebung $w$ in [mm]' , fontsize = 18 )
if len( name[0] ) < 2:
    plt.ylabel( '$\sigma_c[N/mm^2]$' , fontsize = 20 )
    plt.title( 'Zugversuche an stahlfaserbewehrtem Beton mit einer Sollbruchstelle', fontsize = 18 )

else:
    plt.ylabel( '$P$ $[N]$' )
    plt.title( '$Pull-Out$', fontsize = 16 )




Cbsf = CBShortFiber()
w = linspace( 0, 25, 10000 )
all_resp = Cbsf( w, 2.3, 20, 2 * 0.075, 200000, 20, 0, 1, 0, 0, 1e15 )
#plt.plot( w, all_resp, 'k' )
all_resp2 = Cbsf( w, 2.3, 4.5, 2 * 0.075, 200000, 4.5, 0, 1, 0, 0, 1e15 )
#plt.plot( w, all_resp2, 'k' )
#plt.plot( list_x, list_y, color )
#name = ['Pull-Out Test1', 'Pull-Out Test2', 'Fitting']
#plt.legend( name, loc = 'lower right' )
plt.ylim( 0 )
plt.xlim( 0 )
plt.show()



