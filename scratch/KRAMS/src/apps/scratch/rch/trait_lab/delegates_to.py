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
# Created on Sep 11, 2009 by: rch

from enthought.traits.api import \
    HasTraits, Int, DelegatesTo, Instance, PrototypedFrom, on_trait_change, Property

class Worker( HasTraits ):
    
    keep_a_number = Int(14)
    
    @on_trait_change('keep_a_number')
    def _number_changed(self):
        print 'The number has changed'
    
    derived_number = Property( depends_on = 'keep_a_number' )
    def _get_derived_number(self):
        return self.keep_a_number * 10
    
class Chief( HasTraits ):
    
    worker = Instance(Worker)
    
    keep_a_number  = DelegatesTo( 'worker' )

    reuse_a_number = PrototypedFrom( 'worker', prefix = 'keep_a_number' )
    
w = Worker()
c = Chief( worker = w )
print 'retrieved nuber', c.keep_a_number
print 'reused nuber', c.reuse_a_number
c.keep_a_number = 40
print 'retrieved nuber', c.keep_a_number
print 'retrieved nuber', w.keep_a_number

c.reuse_a_number = 42
print 'reused nuber', c.reuse_a_number
print 'reused nuber', w.keep_a_number
    