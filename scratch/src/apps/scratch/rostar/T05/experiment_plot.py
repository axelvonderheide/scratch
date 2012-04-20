'''
Created on Jun 7, 2010

@author: rostislav
'''

from matplotlib import pyplot as p
from numpy import array

def glass():
    # tests with spools
    sp_length = array( [400., 550., 750.] )
    sp_strength = array( [ 766.33, 704.53, 614.67] )
    p.plot( sp_length, sp_strength, 'x--', color = 'black', linewidth = 2,
            label = 'capstan grips' )

    # tests with additional clamps
    clamp_length = array( [0., 50., 100.] )
    clamp_strength = array( [868.47, 951.13 , 1038.03] )
    p.plot( clamp_length, clamp_strength, 's--', color = 'black', linewidth = 2,
            label = 'manual additional clamps' )

    # tests in resin-blocks
    resin_length = array( [100., 200., 300.] )
    resin_strength = array( [1055.92, 913.18 , 928.72 ] )
    p.plot( resin_length, resin_strength, '^--', color = 'black', linewidth = 2,
            label = 'epoxy resin' )

    # tests with statimat
    statimat_std = array( [39.87, 31.48, 64.22, 51.16] )
    statimat_length = array( [135., 300., 500., 800] )
    statimat_strength = array( [795.16, 729.07, 702.93, 677.61] )
    p.plot( statimat_length, statimat_strength, 'o', \
           color = 'blue', linestyle = '--', linewidth = 2,
           label = 'statimat 4U' )
    p.errorbar( elinewidth = 2, x = statimat_length, y = statimat_strength, \
               yerr = statimat_std, color = 'blue', fmt = 'o' )

    cs_area = 0.89 # mm^2
    # tests with statimat
    statimat_4ux_force_std = array( [41.64, 29.63, 37.17, 27.55 ] )
    statimat_4ux_std = statimat_4ux_force_std / cs_area
    # the last one I did not note [rch]
    statimat_4ux_length = array( [50., 150., 300., 500] )
    # force [N]
    statimat_4ux_force = array( [962.61,
                                 904.59,
                                 881,
                                 867] )
    statimat_4ux_strength = statimat_4ux_force / cs_area

    p.plot( statimat_4ux_length, statimat_4ux_strength, 'o', \
           color = 'red', linestyle = '--', linewidth = 3,
           label = 'statimat 4UX' )
    p.errorbar( elinewidth = 2, x = statimat_4ux_length, y = statimat_4ux_strength, \
               yerr = statimat_4ux_std, color = 'red', fmt = 'o' )
    p.ylim( ymin = 0 )
    p.xlabel( 'test length [mm]' )
    p.ylabel( 'breaking stress [MPa]' )
    p.grid()
    p.legend( loc = 'best' )
    p.title( 'AR-glass 2400 tex' )

def carbon():
    # tests with spools
    sp_length = array( [400., 550., 750.] )
    sp_strength = array( [ 1440.56, 1454.78, 1423.45] )
    p.plot( sp_length, sp_strength, 'x--', color = 'black', linewidth = 2 )

    # tests with additional clamps
    clamp_length = array( [0., 50., 100.] )
    clamp_strength = array( [1955.75, 1758.22, 1687.74] )
    p.plot( clamp_length, clamp_strength, 's--', color = 'blue', linewidth = 2 )

    # tests in resin-blocks
    resin_length = array( [100., 200., 300.] )
    resin_strength = array( [1874.99, 1712.21, 1606.77] )
    p.plot( resin_length, resin_strength, '^--', color = 'black', linewidth = 2 )

    # tests with statimat
    statimat_std = array( [65.75, 51.58, 65.68, 65.95] )
    statimat_length = array( [135., 300., 500., 800] )
    statimat_strength = array( [1890.84, 1748.41, 1569.08, 1592.72] )
    statimat_shift = statimat_length - 130

    p.xlabel( 'test length [mm]' )
    p.ylabel( 'breaking stress [MPa]' )
    p.plot( statimat_length, statimat_strength, 'o', \
           color = 'r', linestyle = '--', linewidth = 3 )
    p.plot( statimat_shift, statimat_strength, 'o', \
           color = 'green', linestyle = '--', linewidth = 1 )
    p.errorbar( elinewidth = 2, x = statimat_length, y = statimat_strength, \
               yerr = statimat_std, color = 'r', fmt = 'o' )
    p.errorbar( elinewidth = 1, x = statimat_shift, y = statimat_strength, \
               yerr = statimat_std, color = 'green', fmt = 'o' )
    p.ylim( ymin = 0 )
    p.grid()
    p.legend( ['capstan grips', 'homogenization clamps', 'epoxy resin', 'statimat 4U', 'shift 130 mm'],
              loc = 4 )
    p.title( 'carbon' )

#carbon()
glass()
p.show()
