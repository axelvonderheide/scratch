#-------------------------------------------------------------------------------
#
# Copyright (c) 2009, IMB, RWTH Aachen.
# All rights reserved.
#
# This software is provided without warranty under the terms of the BSD
# license included in simvisage/LICENSE.txt and may be redistributed only
# under the conditions described in the aforementioned license.  The license
# is also available online at http://www.simvisage.com/licenses/BSD.txt
#
# Thanks for using Simvisage open source!
#
# Created on Feb 25, 2010 by: rch

import matplotlib.pyplot as plt
import pickle

from ibvpy.mats.matsXD.matsXD_cmdm.matsXD_cmdm import \
    MATSXDMicroplaneDamage, PhiFnStrainSoftening, PhiFnGeneral, PhiFnStrainHardening 
from matplotlib import rc
    
#phi_fn = PhiFnStrainHardening(Epp = 0.0001, Efp =0.0001, Dfp = 0.4, Elimit = 0.01 )

def plot_tests():
    Epp = 0.0001
    Dfp = 0.4
    Efr = 0.5
    Efp = - ( -1 + Dfp ) * Epp / ( Dfp - Efr ) 
    phi_fn = PhiFnStrainHardening(Epp = Epp, Efp = Efp, Dfp = Dfp, Elimit = 0.006 )
    
    legend = [ '7-a', '9-a', '9-u' ]
    
    for layout in legend:
        
        file_name = '%s_MAG-03-07.mats' % layout
        file = open( file_name, 'r' )
        mfn = pickle.load( file )
        file.close()
    
        plt.plot(mfn.xdata, mfn.ydata)
    
    plt.plot( phi_fn.mfn.xdata, phi_fn.mfn.ydata )
    plt.legend( legend + ['continuous'] )
    plt.title('damage functions')
    plt.ylabel('$\phi$')
    plt.xlabel('$\varepsilon$')
    plt.axes.set_axis_bgcolor(color = 'white')
    plt.show()

def plot_phi_fn():
    
    rc('text', usetex=True)
    rc('axes', labelsize='32' )
    rc('xtick', labelsize = 16 )
    rc('ytick', labelsize = 16 )
    
    phi_fn = PhiFnStrainHardening(Epp = 0.0001, Efp =0.01, Dfp = 0.4, Elimit = 0.006 )

    plt.plot( phi_fn.mfn.xdata, phi_fn.mfn.ydata, color = 'black', linewidth = 2 )
    plt.title('damage function')
    plt.ylabel('$\phi$ [-]')
    plt.xlabel('$\epsilon$ [-]')
    plt.grid( True )
    #plt.set_axis_bgcolor(color = 'white')
    plt.show()
plot_phi_fn()
