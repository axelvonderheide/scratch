'''target curve for the fitting procedure of the phi-function:
linear elastic response for the plane stress state (E = 34000)
(= body can be deformed without affecting the stresses in the 1-direction)
for plane_strain use E/(1-nu**2) instead!:
'''
from numpy import array

# first principle strains:
xdata = array( [0., 1.0] )

# first principle stresses:
ydata = array( [0., 34000] )


