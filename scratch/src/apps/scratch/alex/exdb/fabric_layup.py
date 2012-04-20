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
    DelegatesTo, Expression

from fabric_layout import \
    FabricLayout
    
from promod.simdb.simdb_class import \
    SimDBClass, SimDBClassExt

# @todo: it should probably specify the number 
# of layers per meter/centimeter thickness, otherwise, there is 
# an implicit information about the thickness of the cross section.
# 
class FabricLayup( SimDBClass ):
    '''Describes the type and arrangement of the textile reinforcement
    '''    
    # NOTE: In the specification of a textile cross section:
    # - only one type of textile can be used.
    # - the amount of cross sectional area 'a_tex' is captured as a 
    #   constant smeared value which is calculated based on the
    #   percentage of the textile layers which have an orientation 
    #   in 0- or 90-degree direction with the 1-direction of the structure,
    #   'r_tex_0' and 'r_tex_90' respectively.
    
    # technical textile used
    #
    fabric_layout = Instance( FabricLayout )
    
    # Specify the equidistant spacing in thickness direction between the textile layers
    s_tex_z  = Float( unit = 'm'    , input = True, auto_set = False, enter_set = True )

    # percentage of the orientation in 0- and 90-degree, respectively.
    r_tex_0  = Float( unit = '-', input = True, auto_set = False, enter_set = True )
    r_tex_90 = Float( unit = '-', input = True, auto_set = False, enter_set = True )
    
#    def get_n_layers( self, thickness ):
#        return int( thickness / self.s_tex_z )   #[m/m]
#    
#    def get_n_layers_0( self, thickness ):
#        n_layers = self.get_n_layers( thickness )
#        return int( self.r_tex_0 * n_layers )
#    
#    def get_n_layers_90( self, thickness ):
#        n_layers = self.get_n_layers( thickness )
#        return int( self.r_tex_90 * n_layers )

# Setup the database class extension 
#
FabricLayup.db = SimDBClassExt(
            klass = FabricLayup,
            constants = { 
              '7-a_MAG-03-07'       : FabricLayup (
                                           r_tex_0     = 4./7,
                                           r_tex_90    = 3./7,
                                           s_tex_z     = 0.030/(7+1),
                                           fabric_layout = FabricLayout.db['MAG-07-03']
                                           ),
              '7-u_MAG-03-07'       : FabricLayup (
                                           r_tex_0     = 1.,
                                           r_tex_90    = 0.,
                                           s_tex_z     = 0.030/(7+1),
                                           fabric_layout = FabricLayout.db['MAG-07-03']
                                           ),
              '7-u-all90_MAG-03-07' : FabricLayup (
                                           r_tex_0     = 0.,
                                           r_tex_90    = 1.,
                                           s_tex_z     = 0.030/(7+1),
                                           fabric_layout = FabricLayout.db['MAG-07-03']
                                           ),
              '9-a_MAG-03-07'       : FabricLayup (
                                           r_tex_0     = 4./9,
                                           r_tex_90    = 5./9,
                                           s_tex_z     = 0.030/(9+1),
                                           fabric_layout = FabricLayout.db['MAG-07-03']
                                           ),
              '9-u_MAG-03-07'       : FabricLayup (
                                           r_tex_0     = 1.,
                                           r_tex_90    = 0.,
                                           s_tex_z     = 0.030/(9+1),
                                           fabric_layout = FabricLayout.db['MAG-07-03']
                                           ),
              '9-u-all90_MAG-03-07' : FabricLayup (
                                           r_tex_0     = 0.,
                                           r_tex_90    = 1.,
                                           s_tex_z     = 0.030/(9+1),
                                           fabric_layout = FabricLayout.db['MAG-07-03']
                                           ),
              '10-a_2D-02-06a'      : FabricLayup (
                                           r_tex_0     = 0.5,
                                           r_tex_90    = 0.5,
                                           s_tex_z     = 0.030/(10+1),
                                           fabric_layout = FabricLayout.db['2D-02-06a']
                                           ),
              '10-u_2D-02-06a'      : FabricLayup (
                                           r_tex_0     = 1.,
                                           r_tex_90    = 0.,
                                           s_tex_z     = 0.030/(10+1),
                                           fabric_layout = FabricLayout.db['2D-02-06a']
                                           ),
              '8-u_2D-02-06a'       : FabricLayup (
                                           r_tex_0     = 0.5,
                                           r_tex_90    = 0.5,
                                           s_tex_z     = 0.030/(8+1),
                                           fabric_layout = FabricLayout.db['2D-02-06a']
                                           ),
          })

##------------------------------------------------------------------
## layup view:
##------------------------------------------------------------------
#
#layup_view = View(
#                  Item('object.fabric_layout_key',),
#                  Item('object.n_layers'),
#                  Item('object.orientation'),
#                  resizable = True,
#                  scrollable = True
#                  )

if __name__ == '__main__':
    FabricLayup.db.configure_traits()