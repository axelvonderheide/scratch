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
# Created on Jul 23, 2010 by: rch

'''
Created on Jul 23, 2010

@author: jakub
'''

import numpy as np
import matplotlib.pylab as plt

N = 4
#jakub#rch#alex#andy
orderMeans = np.array( [2. * 8.8 + 2 * 10., 4. * 8.8, 10., 4. * 8.8] )
orderStd = ( 2, 3, 4, 1 )

ind = np.arange( N )  # the x locations for the groups
width = 0.3       # the width of the bars

fig = plt.figure()
ax = fig.add_subplot( 111 )
rects1 = ax.bar( ind, orderMeans, width, color = 'r', yerr = orderStd )

creditMeans = np.array( [85.6 + 3. * 3.3, 4 * 3.3, 2. * 3. + 3.3, 0], dtype = float )
creditStd = ( 3., 5., 2., 3. )
rects2 = ax.bar( ind + width, creditMeans, width, color = 'y', yerr = creditStd )

resultMeans = orderMeans - creditMeans

resultStd = orderStd
rects3 = ax.bar( ind + 2 * width, resultMeans, width, color = 'b', yerr = resultStd )

# add some
ax.set_ylabel( r'EURO' )
ax.set_title( 'Bewerage Optimizing Problem' )
ax.set_xticks( ind + width )
ax.set_xticklabels( ( 'Jakub', 'RCH', 'Alex', 'Andy' ) )

ax.legend( ( rects1[0], rects2[0], rects3[0] ), ( 'Order', 'Credit', 'Result' ) )

def autolabel( rects, values ):
    # attach some text labels
    for  rect, val in zip( rects, values ):
        height = rect.get_height()
        if val > 0:
            ax.text( rect.get_x() + rect.get_width() / 2., height + 0.4, '%0.1f' % val,
                    ha = 'center', va = 'bottom' )
        else:
            ax.text( rect.get_x() + rect.get_width() / 2., 0.4, '%0.1f' % val,
                    ha = 'center', va = 'bottom' )
autolabel( rects1, orderMeans )
autolabel( rects2, creditMeans )
autolabel( rects3, resultMeans )
plt.show()
