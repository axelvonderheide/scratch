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
# Created on Mar 18, 2011 by: rch

from numpy import loadtxt, argmax, array, arange, zeros, average, min, max

from promod.simdb import SimDB

from os.path import join

import pylab as p

def plot_csv_arr( ax, file, ix_arr, p, delimiter = ',',
                  xfactor = 1.0, yfactor = 1.0,
                  color = 'black' ):
    '''
    '''
    converters = {}
    for ix in ix_arr.flatten():
        converters[ix] = s_or_0

    ep_arr = loadtxt( file, delimiter = delimiter,
                      dtype = float,
                      skiprows = 21,
                      usecols = ix_arr.flatten(),
                      converters = converters )

    # get the number of columns
    #
    n_curves = ep_arr.shape[1] / 2
    ix_arr = arange( ep_arr.shape[1] ).reshape( n_curves, 2 )

    linestyles = ['dashed', 'solid', 'dotted' ]

    ymax_list = []
    for i, j in ix_arr:
        idx = j / 2
        max_idx = argmax( ep_arr[:, j] )
        ydata = yfactor * ep_arr[:max_idx, j]
        xdata = xfactor * ep_arr[:max_idx, i]
        ax.plot( xdata, ydata , color = color,
                 linestyle = linestyles[idx],
                 label = 'Versuch %d' % ( idx + 1 ) )
        ymax_list.append( ydata[-1] )

    ax.legend( loc = 'lower right' )
    ymax_arr = array( ymax_list, dtype = float )

    ymin, ymax = min( ymax_arr ), max( ymax_arr )

    ax.axhspan( ymin, ymax, facecolor = '0.5', alpha = 0.5 )

    xdelta = ax.get_xlim()[-1]
    print xdelta
    ydelta = 30
    ax.text( 0.1 * xdelta, ydelta + ymax, 'max = %4.0f N' % ymax )
    ax.text( 0.1 * xdelta, -ydelta + ymin, 'min = %4.0f N' % ymin, va = 'top' )

    return average( ymax_arr ), ymin, ymax

if __name__ == '__main__':

    from promod.simdb import SimDB
    simdb = SimDB()

    s_or_0 = lambda s: float( s or 0 )

    # Plot the pullout tests
    #

    p.subplots_adjust( wspace = 0.1 )

    ax1 = p.subplot( 1, 2, 2 )
    ax1.set_xlabel( 'Weg [mm]' )
    data_dir = join( simdb.exdata_dir, 'trc', 'yarn_pull_out', 'ep_ypo' )

    ep_ytt_file = join( data_dir, 'Pull-Out-AR epstfstd.csv' )
    ix_arr = array( [[0, 1], [5, 6], [10, 11]] )
    ypo_mean, ypo_min, ypo_max = plot_csv_arr( ax1, ep_ytt_file, ix_arr, p )

    ax2 = p.subplot( 1, 2, 1, sharey = ax1 )
    ax2.set_xlabel( 'Dehnung [-]' )
    ax2.set_ylabel( 'Kraft [N]' )

    # Plot the tensile tests
    #
    data_dir = join( simdb.exdata_dir, 'trc', 'yarn_tensile_tests', 'ep_ytt' )

    ep_ytt_file = join( data_dir, 'ARG GZ Epoxies 125mm.csv' )
    ix_arr = array( [[0, 1], [3, 4], [6, 7]] )
    ytt_125_mean, ytt_125_min, ytt_125_max = plot_csv_arr( ax2, ep_ytt_file, ix_arr, p, yfactor = 0.89 )

#    ep_ytt_file = join( data_dir, 'ARG GZ Epoxies 500mm.csv' )
#    ix_arr = array( [[0, 1], [3, 4], [6, 7], [9, 10]] )
#    ytt_500_mean, ytt_500_min, ytt_500_max = plot_csv_arr( ep_ytt_file, ix_arr, p, yfactor = 0.89 )

#    print 1 - ypo_mean / ytt_125_mean
#    print ypo_mean, ytt_125_mean, ytt_500_mean
#    print ypo_min, ytt_125_min, ytt_500_min
#    print ypo_max, ytt_125_max, ytt_500_max

    ax1.set_title( 'Garn-Pull-Out' )
    ax1.grid( True )
    ax2.set_title( 'Garn-Zugversuch' )
    ax2.grid( True )
    p.setp( ax1.get_yticklabels(), visible = False )

    p.ylim( 0, 2000 )
    p.show()
