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
# Created on Mar 11, 2011 by: rch

from sympy import *

x = Symbol( 'x' )
y = Symbol( 'y' )
print solve( x ** 2 - y, x )

var( 'beta k la' )

expr = k / ( 1 + k ) * beta - la

print expr

k_solve = solve( expr, [k] )

print k_solve

Plot( k_solve[0], [la, 0, 10], [beta, 0, 10] )
