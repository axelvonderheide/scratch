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
# Created on Sep 9, 2009 by: rch
from math import fabs

def get_one_if_same_sign( a, b ):
    sa = fabs(a) / a
    sb = fabs(b) / b
    return fabs( 1./2. * (sa + sb) )

print 'yes', get_same_sign( 1., 4. )
print 'no', get_same_sign( 1., -4. )
print 'no', get_same_sign( -1., 4. )
print 'yes', get_same_sign( -1., -4. )