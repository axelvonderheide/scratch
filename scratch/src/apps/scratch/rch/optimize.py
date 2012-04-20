
from scipy.optimize import fsolve

from scipy import stats

from numpy import array

from math import sqrt
#    distr_type = Enum('norm','gamma','weibull_min','uniform')


def get_distr_residuum( params, mean, stdev, skew ):
    print 'maan',mean,'variance',stdev
    print params
    shape = params[0]
    loc   = params[1]
    scale = params[2]
    mom_eq = stats.weibull_min( shape, loc = loc, scale = scale, moments = 'mvs' )
    moms = mom_eq.stats()
    print '---------------------------------'
    print moms
    mean_resid     = moms[0] - mean
    stdev_resid = sqrt( moms[1] ) - stdev
    skew_resid     = moms[2] - skew
    print 'mean', moms[0], 'stdev', sqrt( moms[1] )
    return ( mean_resid, stdev_resid, skew_resid )
    #return array( [ mean_resid, stdev_resid, skew_resid ], dtype = 'float_' )

#print fn( [2, 0, 1], 0.0, 1.0 )

res = fsolve( get_distr_residuum, [1,1,1], 
              args = (100., 10.0, 1.0 ) )
print 'xxxx'
print res

#print fn( [5,0,10], 0, 0, 0 )
