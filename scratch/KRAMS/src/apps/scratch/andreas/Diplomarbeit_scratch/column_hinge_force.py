.0726
'''
Created on Jul 30, 2010

@author: abach
'''
from enthought.traits.api import \
    HasTraits, Directory, List, Int, Float, Any, Enum, \
    on_trait_change, File, Constant, Instance, Trait, \
    Array, Str, Property, cached_property, WeakRef, \
    Dict, Button, Color, Bool

from scipy import integrate
from numpy import *
from matplotlib import pyplot
from math import pi

#steps = int ( 1000 )
#x_i = linspace( 0, 3, steps ) # m
#M_0 = linspace( 72,6, 94,9, steps )  #m
#M_1 = linspace(0,1,steps)
#I_i = ( ( -0.1 * x_i / 3.0 + 0.45 ) ** 4 ) / 12  # m^4
#E = 38800  # MN/m^2
#
#Func = M ** 2 / ( E * I_i ) # 1/ m
#Integral = integrate.trapz( Func, x_i )
#Force = 0.0052 / Integral

def func_u_z(width_bottom, width_top):
    steps = int ( 1000 )
    x_i = linspace( 0, 3, steps ) # m
    M_0 = linspace( 1.0,1.0, steps )  #m
#    M_0 = linspace( .0726+0.02688,0.0949+0.02688, steps )  #m

    M_1 = linspace( 1.0,1.0, steps )
    I_i = ( (  width_top + x_i * (width_bottom - width_top)/3 ) ** 4 ) / 12  # m^4
    E = 38800  # MN/m^2
    
    Func = M_0 * M_1 / ( E * I_i ) # 1/ m
    Integral = integrate.trapz( Func, x_i )
    
    return Integral *4

def func_d_u(width_bottom, width_top):
    steps = int ( 1000 )
    x_i = linspace( 0, 3, steps ) # m
    N_0 = linspace( 0.040,0.040, steps )  #m
#    M_0 = linspace( .0726+0.02688,0.0949+0.02688, steps )  #m

    N_1 = linspace( 1.0,1.0, steps )
    A_i = ( (  width_top + x_i * (width_bottom - width_top)/3 ) ** 2 )   # m^4
    E = 38800  # MN/m^2
    
    Func = N_0 * N_1 / ( E * A_i ) # 1/ m
    Integral = integrate.trapz( Func, x_i )
    
    return Integral *4


width_bottom = arange(0.35,0.901,0.05)

print func_u_z(0.45,0.35)*26.4
print func_d_u(0.45,0.35)*1000
#u_z_0=[]
#for i in width_bottom:
#    u_z_0 = u_z_0 + [func_u_z(i,0.45)]
#
#u_z_1=[]
#for i in width_bottom:
#    u_z_1 = u_z_1 + [func_u_z(i,0.50)]
#    
#u_z_2=[]
#for i in width_bottom:
#    u_z_2 = u_z_2 + [func_u_z(i,0.55)]
#
#fig = pyplot.figure( facecolor = "white" )
#ax1=fig.add_subplot( 1, 1, 1 )
#ax1.plot( width_bottom , u_z_0 , color ='blue', label='Breite Kopf .45')
#ax1.plot( width_bottom , u_z_1 , color ='red', label='Breite Kopf .50')
#ax1.plot( width_bottom , u_z_2 , color ='green', label='Breite Kopf .55')
#ax1.set_ylabel( 'u_z [m]', fontsize = 22)
#ax1.set_xlabel( 'Breite Stuetzenfuss [m] ', fontsize = 22)
#pyplot.legend()
#pyplot.show()

#fig.savefig( filename,orientation = 'portrait', bbox_inches='tight' )
#pyplot.clf()