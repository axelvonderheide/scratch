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
# Created on Feb 23, 2010 by: rch

from enthought.traits.api import \
    HasTraits, Directory, List, Int, Float, Any, \
    on_trait_change, File, Constant, Instance, Trait, \
    Array, Str, Property, cached_property, WeakRef, \
    Dict, Button, Bool, Enum, Event, implements, \
    DelegatesTo

from promod.simdb.simdb_class import \
    SimDBClass, SimDBClassExt
    
class FabricLayout( SimDBClass ):
    '''Comprises the characteristics of the textile reinforcement
    '''
    # technical name of the reinforcement:
    reinf_name = Str

    # cross sectional area of the reinforcement [mm^2/m]: 
    # in 0 and 90 degree orientation: 
    a_tex_0  = Float(unit = 'mm^2/m', simdb = True )
    a_tex_90 = Float(unit = 'mm^2/m', simdb = True )
    E_tex_0  = Float( unit = 'MPa' ,  simdb = True )
    E_tex_90 = Float( unit = 'MPa' ,  simdb = True )

    # spacing of the textile mesh [mm]: 
    s_tex_0  = Float(unit = 'mm', simdb = True )
    s_tex_90 = Float(unit = 'mm', simdb = True )

#    # cross sectional area of one roving in 0-direction [mm^2]: 
#    A_roving_0  = Property( Float, unit = 'mm^2', simdb = True  )
#    def _get_A_roving_0( self ):
#        return self.a_tex_0  * self.s_tex_0/1000
#
#    # cross sectional area of one roving in 90-direction [mm^2]: 
#    A_roving_90 = Property( Float, unit = 'mm^2', simdb = True  )
#    def _get_A_roving_90( self ):
#        return self.a_tex_90 * self.s_tex_90/1000 
    
#    # Discrete number of rovings used in the specimen in 
#    # longitudinal direction for one layer of textile fabric of size 'width'
#    # distiguished in 0- and  90-direction, respectively. 
#    # The value is rounded to the next smaller integer. 
#    def get_n_rovings_0( self, width ):
#        width = width * 1000.               #[mm] 
#        return int( width / self.s_tex_0 )  #[mm/mm]
#
#    def get_n_rovings_90( self, width ):
#        width = width * 1000.               #[mm]       
#        return int( width / self.s_tex_90 ) #[mm/mm]

# Setup the database class extension 
#
FabricLayout.db = SimDBClassExt(
            klass = FabricLayout,
            constants = {         
               'MAG-07-03' : FabricLayout(
                                           reinf_name = 'MAG-07-03',
                                           a_tex_0  = 107.89,
                                           a_tex_90 = 106.61,
                                           E_tex_0  = 70000,
                                           E_tex_90 = 70000,
                                           s_tex_0  = 8.3,
                                           s_tex_90 = 8.4,
                                           ),
               '2D-02-06a' : FabricLayout(
                                           reinf_name = '2D-06-06a',
                                           a_tex_0  = 71.65,
                                           a_tex_90 = 53.31,
                                           E_tex_0  = 70000,
                                           E_tex_90 = 70000,
                                           s_tex_0  = 12.5,
                                           s_tex_90 = 8.4,
                                           )
             }
            )

if __name__ == '__main__':
    FabricLayout.db.configure_traits()