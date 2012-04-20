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
# Created on Jan 29, 2010 by: rch

from math import sqrt
from enthought.traits.api import HasTraits, Float, Property

class Coord( HasTraits ):
    
    x = Float( 0.0 )
    y = Float( 0.0 )
    
    length = Property( Float, depends_on = 'x,y' )
    def _get_length(self):
        return sqrt( self.x**2 + self.y**2 )
    
    def __init__(self, **kw):
        print kw
        super( HasTraits, self ).__init__(**kw)
    
c = Coord( x = 1.0, y = 2.0 )
c.configure_traits()

