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
# Created on Jun 30, 2010 by: rch

import numpy as np
import matplotlib.pyplot as plt

fig = plt.figure()
ax1 = fig.add_subplot( 111 )
t = np.arange( 0.01, 10.0, 0.01 )
s1 = np.exp( t )
ax1.plot( t, s1, 'b-' )
ax1.set_xlabel( 'time (s)' )
# Make the y-axis label and tick labels match the line color.
ax1.set_ylabel( 'exp', color = 'b' )
for tl in ax1.get_yticklabels():
    tl.set_color( 'b' )


ax2 = ax1.twinx()
s2 = np.sin( 2 * np.pi * t )
ax2.plot( t, s2, 'r.' )
ax2.set_ylabel( 'sin', color = 'r' )
for tl in ax2.get_yticklabels():
    tl.set_color( 'r' )
plt.show()
