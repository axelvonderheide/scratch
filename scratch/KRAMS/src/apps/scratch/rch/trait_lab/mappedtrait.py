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
# Created on Aug 24, 2009 by: rch

from enthought.traits.api import HasTraits, Callable, Trait, Enum, on_trait_change, Array, Float, Str

class MyModel( HasTraits ):

    name = Str
    input_data = Array( float )

    def calculate( self ):
        print 'do something'


mymodel = MyModel()
mymodel.configure_traits()
