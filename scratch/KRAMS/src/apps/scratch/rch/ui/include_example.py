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

from enthought.traits.api import HasTraits, Str, Int
from enthought.traits.ui.api import View, Include, Group

class Person(HasTraits):
    name = Str
    age = Int
    pview = Group( 'name', 'age' )
    person_view = View('name', Include('extra'), 'age', kind='modal')

class LocatedPerson(Person):
    street = Str
    city = Str
    state = Str
    zip = Int
    extra = Group('street', 'city', 'state', 'zip')
    moreview = View( Include('extra'), Include('person_view') )

lp = LocatedPerson()
lp.configure_traits( view = 'moreview' )
