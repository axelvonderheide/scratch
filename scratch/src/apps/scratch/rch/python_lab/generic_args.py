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
# Created on Jun 17, 2010 by: rch

# How to pass extensible set of arguments through three levels of functions

from enthought.traits.api import HasTraits, Callable

def mats_f1( eps ):
    print eps

def mats_f2( eps, eps_avg ):
    print eps, eps_avg

class fets( HasTraits ):

    mats = Callable

    def __call__( self, eps, *args, **kw ):
        self.mats( eps, *args, **kw )

if __name__ == '__main__':

    # without averaging
    f = fets()
    f.mats = mats_f1
    f( 0.5 )

    # with averaging
    f.mats = mats_f2
    f( 0.5, 0.7 )

    f_dict = {'eps_avg' : 0.9 }
    f( 0.5, **f_dict )



