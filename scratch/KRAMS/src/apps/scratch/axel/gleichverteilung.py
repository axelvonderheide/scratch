'''
Created on 25.09.2011

@author: axel
'''
import numpy as np
from matplotlib import pyplot as plt

x = np.linspace( 0, 10, 100 )
y = 100 - x * 10
plt.plot( x, y , 'k' )
plt.title( 'Prozentualer Anteil an Fasern mit minimaler Einbindetiefe von x', fontsize = 14 )
plt.xlabel( 'Abstand zum Risspunkt in [mm]' , fontsize = 13 )
plt.ylabel( 'Anteil an Fasern in [%]', fontsize = 13 )
plt.show()
