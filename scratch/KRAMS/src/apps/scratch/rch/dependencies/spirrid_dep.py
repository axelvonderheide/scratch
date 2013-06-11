#-------------------------------------------------------------------------------
#
# Copyright (c) 2009, IMB, RWTH Aachen.
# All rights reserved.
#
# This softwe_arare is provided without warranty under the terms of the BSD
# license included in simvisage/LICENSE.txt and may be redistributed only
# under the conditions described in the aforementioned license.  The license
# is also available online at http://www.simvisage.com/licenses/BSD.txt
#
# Thanks for using Simvisage open source!
#
# Created on Sep 8, 2011 by: rch

from enthought.traits.api import \
    HasTraits, Array, Property, Float, cached_property, Callable, Tuple, \
    List, Str, Enum, Int, Instance, Trait, WeakRef, Bool, on_trait_change

from stats.spirrid.code_gen import \
    CodeGenNumpy, CodeGenC
from stats.spirrid.sampling import \
    FunctionRandomization, TGrid, PGrid, MonteCarlo, LatinHypercubeSampling, orthogonalize

import numpy as np

import platform
if platform.system() == 'Linux':
    from time import time as sysclock
elif platform.system() == 'Windows':
    from time import clock as sysclock

#===============================================================================
# Generic implementation of the integral
#===============================================================================
class SPIRRID(FunctionRandomization):
    '''Set of parallel independent responses with random identical distributions.
    '''

    #===========================================================================
    # type of the sampling of the random domain
    #===========================================================================
    sampling_type = Trait('T-grid', {'T-grid' : TGrid,
                                     'P-grid' : PGrid,
                                      },
                          input_change = True)

    sampling = Property(depends_on = '+input_change')
    @cached_property
    def _get_sampling(self):
        print '... recalculating result ...',
        return self.sampling_type_(randomization = self)

