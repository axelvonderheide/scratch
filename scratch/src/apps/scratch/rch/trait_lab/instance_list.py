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
# Created on Aug 7, 2009 by: rchx

# instance_editor_selection.py -- Example of an instance editor
#                                 with instance selection

from enthought.traits.api    \
    import HasStrictTraits, Int, Instance, List, Regex, Str
from enthought.traits.ui.api \
    import View, Item, InstanceEditor
from enthought.traits.ui.instance_choice import \
    InstanceFactoryChoice
from either_type import EitherType

class Person ( HasStrictTraits ):
    name  = Str
    age   = Int
    phone = Regex( value = '000-0000',
                   regex = '\d\d\d[-]\d\d\d\d' )

    traits_view = View( 'name', 'age', 'phone', resizable = True )

class Pet ( HasStrictTraits ):
    name  = Str
    age   = Int

    traits_view = View( 'name', 'age', resizable = True )


class Team ( HasStrictTraits ):

    name    = Str
    captain = EitherType( klasses = [Person, Pet],
                          names = ['person','zviratko'] )

    traits_view = View( Item('name'),
                        Item('_'),
                        Item( 'captain',
                              label='Team Captain',
                              style = 'custom',
                             ),
                        resizable = True,
                        buttons = ['OK'])

if __name__ == '__main__':
    Team( name    = 'Vultures').configure_traits()