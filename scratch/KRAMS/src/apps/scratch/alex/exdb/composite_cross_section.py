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
    DelegatesTo, Callable

from enthought.traits.ui.api import \
    View, Item, VGroup, Group, Spring, HGroup

# overload the 'get_label' method from 'Item' to display units in the label
from traits.ui.item import \
    Item
        
from concrete_mixture \
    import ConcreteMixture

from fabric_layout \
    import FabricLayout

from fabric_layup \
    import FabricLayup
        
from promod.simdb.simdb_class import \
    SimDBClass, SimDBClassExt

 
class CompositeCrossSection( SimDBClass ):
    '''Describes the combination of the concrete mixture and textile cross section.
    '''    

    #--------------------------------------------------------------------------------
    # select the concrete mixture from the concrete database:
    #--------------------------------------------------------------------------------

    concrete_mixture_key = Enum( ConcreteMixture.db.keys(), 
                                 simdb = True, input = True, auto_set = False, enter_set = True )
#    concrete_mixture_key = Property( Str, simdb = True )
    
    concrete_mixture_ref = Property( Instance( SimDBClass ) )
    def _get_concrete_mixture_ref(self):
        return ConcreteMixture.db[ self.concrete_mixture_key ]

    #--------------------------------------------------------------------------------
    # select the textile cross section from the textile database:
    #--------------------------------------------------------------------------------

    fabric_layup_key = Enum( FabricLayup.db.keys(), 
                                      simdb = True, input = True, auto_set = False, enter_set = True  )  
    
    fabric_layup_ref = Property( Instance( SimDBClass ) )
    def _get_fabric_layup_ref(self ):
        return FabricLayup.db[ self.fabric_layup_key ]

    fabric_layout = Property( Instance(FabricLayout) )
    def _get_fabric_layout( self ):
        return self.fabric_layup_ref.fabric_layout

    #--------------------------------------------------------------------------------
    # define the cross section layup:
    #--------------------------------------------------------------------------------

    fabric_layup_list = List( Instance( FabricLayup )) 

    fabric_layup_list_keys = List( Str )
    def _get_fabric_layup_list_keys(self ):
        return [ flu.key for flu in self.fabric_layup_list ]
    
    fabric_layup_list_refs = Property( List( Instance( SimDBClass ) ) )
    def _get_fabric_layup_list_refs(self ):
        return [ FabricLayup.db[ flu_key ] for flu_key in self.fabric_layup_list_keys ]

        





    #--------------------------------------------------------------------------------
    # calculated material properties 
    #--------------------------------------------------------------------------------
    
    # NOTE: up to now only a constant reinforcement ratio for each 
    # textile cross section can be specified. The evaluation is prepared based on the
    # spacing of the textile layers and the composite dimensions. The discrete number
    # of layers in z-direction and the discrete number of rovings in y-direction are
    # calculated automatically and rounded to the next lowest possible integer.  
    
    # @todo: do the unreinforced boundary zones need to be subtracted in the calculation of the plate/tensile test?
    

    def get_A_tex( self, thickness, width ):
        '''total cross-sectional-area of the textile reinforcement.
        ''' 
        A_roving_0   = self.fabric_layout.A_roving_0   #[mm^2]    
        A_roving_90  = self.fabric_layout.A_roving_90  #[mm^2]    
        
        n_layers_0   = self.fabric_layup_ref.get_n_layers_0( thickness )
        n_layers_90  = self.fabric_layup_ref.get_n_layers_90( thickness )

        n_rovings_0  = self.fabric_layout.get_n_rovings_0( width )
        n_rovings_90 = self.fabric_layout.get_n_rovings_90( width )

        A_tex = A_roving_0  * n_rovings_0  * n_layers_0 + \
                A_roving_90 * n_rovings_90 * n_layers_90
        
        return A_tex

    def get_rho_c( self, thickness, width ):
        '''Return a function for the reinforcement ratio of the composite 
        material depending on the textile fabric layout and the dimensions
        of the composite cross section in x-direction.
        ''' 
        A_c   = width * thickness * 1000000         #[mm^2]
        A_tex = self.get_A_tex( thickness, width )  #[mm^2]
        rho_c   = A_tex / A_c  #[-]
        return rho_c 

    def get_E_c_time( self, age, thickness, width ):
        '''function for the composite E-modulus depending of the concrete age:
        '''
        E_m = self.concrete_mixture_ref.get_E_m_time( age )
        rho = self.get_rho_c( thickness, width )
        E_tex = self.fabric_layup_ref.fabric_layout.E_tex_0
        return (1-rho) * E_m + rho * E_tex  

    rho_cc = Property( Float, unit = '-' )
    def _get_rho_cc( self ):
        '''Return the reinforcement ratio of the composite 
        of a periodic cell with a given textile cross 
        section in x-direction as a reference value.
        ''' 
        a_tex_0  = self.fabric_layup_ref.fabric_layout.a_tex_0   #[mm^2/mm]
        a_tex_90 = self.fabric_layup_ref.fabric_layout.a_tex_90  #[mm^2/mm]
        s_tex_z  = self.fabric_layup_ref.s_tex_z * 1000.                 #[mm]
        r_tex_0  = self.fabric_layup_ref.r_tex_0                         #[-]
        r_tex_90 = self.fabric_layup_ref.r_tex_90                        #[-]
        a_tex  = a_tex_0 * r_tex_0 + a_tex_90 * r_tex_90                          #[mm^2/mm]
        rho_cc = a_tex / ( s_tex_z * 1000.)                                       #[-]
        return rho_cc 

    E_cc = Property( Float, unit = 'MPa' )
    def _get_E_cc( self ):
        '''function for the composite E-modulus depending of the concrete age:
        '''
        E_m = self.concrete_mixture_ref.E_m28
        rho = self.rho_cc
        E_tex = self.fabric_layup_ref.fabric_layout.E_tex_0
        return (1-rho) * E_m + rho * E_tex  



    #--------------------------------------------------------------------------------
    # view:
    #--------------------------------------------------------------------------------

    traits_view = View( Item('key', style = 'readonly', show_label = False ),
                        VGroup(
                        Group(
                        HGroup( Spring(),
                                Item('concrete_mixture_key', 
                                     style = 'readonly', show_label = False, emphasized = True ),
                                Spring()
                              ),
                        Item('concrete_mixture_ref@', show_label = False ),
                        label = 'concrete mixture'
                        ),
                        Group(
                        HGroup( Spring(),
                                Item('fabric_layup_key', 
                                     style = 'readonly', show_label = False, emphasized = True),
                                Spring()),
                        Item('fabric_layup_ref@', show_label = False ),
                        label = 'fabric layup',
                        ),
                        Group(
                        Item('rho_cc', 
                             style = 'readonly', show_label = True ),
                        Item('E_cc', 
                             style = 'readonly', show_label = True ),
                        label = 'composite cross section',
                        ),

                        #layout = 'tabbed',
                        orientation = 'vertical',
                        ),
                        scrollable = True
                        )

#  Group( 
#                            Item('rho', style = 'readonly', show_label = False ),
#                            Item('rho', style = 'readonly', show_label = False ),
#                            label = 'composite cross section')
#                             ),










# Setup the database class extension 
#
CompositeCrossSection.db = SimDBClassExt(
            klass = CompositeCrossSection,
            constants = { 
              '7-a_MAG-03-07_PZ-0708-1' : CompositeCrossSection (
                                            fabric_layup_key = '7-a_MAG-03-07', 
                                            concrete_mixture_key = 'PZ-0708-1'
                                            ) 
          })

if __name__ == '__main__':
    CompositeCrossSection.db.configure_traits()