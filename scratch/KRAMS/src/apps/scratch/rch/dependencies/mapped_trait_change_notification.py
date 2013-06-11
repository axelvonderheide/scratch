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
# Created on Oct 5, 2011 by: rch

from enthought.traits.api import HasStrictTraits, Trait, on_trait_change, \
    Property, cached_property

from stats.spirrid.sampling import PGrid, TGrid, FunctionRandomization
from stats.spirrid import SPIRRID

class C(FunctionRandomization):

    sampling_type = Trait('T-grid', {'T-grid' : TGrid,
                            'P-grid' : PGrid}, input_change = True)

    sampling = Property(depends_on = '+input_change')
    @cached_property
    def _get_sampling(self):
        print '... recalculating result ...',
        return self.sampling_type_(randomization = self)

c = C()
print 'got', c.sampling
c.sampling_type = 'P-grid'
print 'got', c.sampling
c.sampling_type = 'T-grid'
print 'got', c.sampling
print 'got', c.sampling

print '====='
s = SPIRRID()
print 'got', s.sampling
s.sampling_type = 'P-grid'
print 'got', s.sampling
s.sampling_type = 'T-grid'
print 'got', s.sampling
print 'got', s.sampling
