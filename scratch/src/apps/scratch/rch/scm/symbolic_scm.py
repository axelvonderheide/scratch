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


from sympy import var, Plot, exp

def dataset( **kw ):
    return kw

var( 'sigma_c sigma_Rc m' )

weibull = 1 - exp( -( sigma_c / sigma_Rc ) ** m )

data = dataset( sigma_Rc = 0.5, m = 3.0 )

weibull_data = weibull.subs( data )

Plot( weibull_data, [0, 5] )
