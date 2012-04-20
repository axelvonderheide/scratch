'''
Created on 08.09.2011

@author: axel
'''
from math import pi



bewehrungstahl_durchmesser = 8.
streckgrenze_stahl = 500.


############### Reinforcement ####################################
Er = 200000
Dr = bewehrungstahl_durchmesser
Ar = ( Dr / 2 ) ** 2 * pi
max_sigma_r = streckgrenze_stahl



###########fibers####################
r = 0.075
A_fiber = pi * r ** 2
fiber_length = 9.

############specimen######################################
length = 500
width = 30
height = 30
A_spec = width * height
fiber_volume_fraction = 3
Em = 30000
A_c = A_spec - Ar


############specimen tests#################
crack_sigma_c = 6
A_test = 27 * 27

#tau 
test_crack_dist = 50
A_test = 30 * 30 - 3 ** 2 * pi
max_sigma_concrete = 2

T_r = max_sigma_concrete * A_test / 50
tau = T_r / ( pi * 6 )




#################### Computation ##################################

####estimated E-modulus of composite 

E_c = ( A_c * fiber_volume_fraction / 100 * Er + Em * ( A_c - ( 1 - fiber_volume_fraction / 100 ) ) ) / A_c

####breaking strain of composite 

max_strain_c = crack_sigma_c / E_c



max_eps_r = max_sigma_r / Er
Tr = tau * Dr * pi
crack_dist = ( Er * Ar + E_c * A_c ) * max_strain_c / Tr

print 'Bruchdehnung composite' , max_strain_c * 1000, 'promille'
print 'Bruchdehnung reinf.' , max_eps_r * 1000, 'promille'
print 'Risslast im Pruefquerschnitt' , ( Er * Ar + E_c * A_c ) * max_strain_c, 'N'
print 'Kraft bis zur Streckgrenze' , max_eps_r * Ar * Er, 'N'
print 'Schubfluss pro mm', Tr, 'N'
print 'Rissabstand', crack_dist, 'mm'
print 'voraussichtliche Anzahl an Rissen', length / crack_dist

