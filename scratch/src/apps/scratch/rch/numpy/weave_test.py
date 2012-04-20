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
# Created on Apr 13, 2010 by: rch

from scipy import weave

code = r"""
int i;
// py::tuple results(2);
//for (i=0; i<a.length(); i++) {
//     a[i] = i;
//}
//results[0] = 3.0;
//results[1] = 4.0;
return_val = 10.0 * a;
"""
a = 20
res = weave.inline(code,['a'])

print res