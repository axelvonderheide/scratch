'''
Created on May 11, 2010

@author: rostislav
'''

from math import pi

# work of one screw rotation times torque (torque * 2 * pi) equals
# work of axial screw displacement times axial force (F * P)
# zero friction assumed

torque = 0.25 #Nm
P = 0.001 #lead of screw thread = 1 mm for M6
N = 4 # No of screws
F = N * torque * 2 * pi / P

print 'force of', N, 'screw(s) =', F, 'N'

# stress in MPa

A = 0.05*0.025 # area of the polymer layer
print 'stress =', F/A/10**6, 'MPa'