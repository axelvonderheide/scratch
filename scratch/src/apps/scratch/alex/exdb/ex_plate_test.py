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
# Created on Feb 15, 2010 by: rch

from enthought.traits.api import \
    HasTraits, Directory, List, Int, Float, Any, \
    on_trait_change, File, Constant, Instance, Trait, \
    Array, Str, Property, cached_property, WeakRef, \
    Dict, Button, Bool, Enum, Event, implements, DelegatesTo

from enthought.util.home_directory import \
    get_home_directory

from enthought.traits.ui.api import \
    View, Item, DirectoryEditor, TabularEditor, HSplit, VGroup, \
    TableEditor, EnumEditor, Handler, FileEditor, VSplit, Group, \
    HGroup, Spring
    
# overload the 'get_label' method from 'Item' to display units in the label
from traits.ui.item import \
    Item
    
from enthought.traits.ui.table_column import \
    ObjectColumn
    
from enthought.traits.ui.menu import \
    OKButton, CancelButton
    
from enthought.traits.ui.tabular_adapter \
    import TabularAdapter

from traits.editors.mpl_figure_editor import MPLFigureEditor
from matplotlib.figure import Figure
    
import os

import csv

from numpy import array, fabs, where, copy, ones

from scipy.io import read_array
from numpy import \
    loadtxt, argmax, polyfit, poly1d, frompyfunc, dot

from enthought.traits.ui.table_filter \
    import EvalFilterTemplate, MenuFilterTemplate, RuleFilterTemplate, \
           EvalTableFilter

from enthought.traits.ui.api import TextEditor 
float_editor_9g = TextEditor( evaluate=float, format_str="%.9g", enter_set=True, auto_set=False )
float_editor_0f = TextEditor( evaluate=float, format_str="%.0f", enter_set=True, auto_set=False )
float_editor_1f = TextEditor( evaluate=float, format_str="%.1f", enter_set=True, auto_set=False )
float_editor_2f = TextEditor( evaluate=float, format_str="%.2f", enter_set=True, auto_set=False )
float_editor_3f = TextEditor( evaluate=float, format_str="%.3f", enter_set=True, auto_set=False )
float_editor_4f = TextEditor( evaluate=float, format_str="%.4f", enter_set=True, auto_set=False )

#-- Tabular Adapter Definition -------------------------------------------------

from string import replace
from os.path import exists

#-----------------------------------------------------------------------------------
# ExDesignReader
#-----------------------------------------------------------------------------------
from enthought.traits.ui.file_dialog  \
    import open_file, FileInfo, TextInfo, ImageInfo

from enthought.traits.ui.api \
    import View, Item, TabularEditor
    
from enthought.traits.ui.tabular_adapter \
    import TabularAdapter

from ex_type import ExType

data_file_editor = FileEditor( filter = ['*.DAT'] )

from ex_type import ExType
from i_ex_type import IExType

from mathkit.array.smoothing import smooth

from promod.matdb.trc.fabric_layup \
    import FabricLayup
    
from promod.matdb.trc.fabric_layout \
    import FabricLayout 
  
from promod.matdb.trc.concrete_mixture \
    import ConcreteMixture

from promod.matdb.trc.composite_cross_section import \
    CompositeCrossSection



class ExPlateTest( ExType ):
    '''Read the data from the directory
    '''
    
    implements( IExType )

    #--------------------------------------------------------------------------------
    # specify inputs:
    #--------------------------------------------------------------------------------
    
    edge_length   = Float( 1.25, unit = 'm', input = True, auto_set = False, enter_set = True, table_field = True )
    thickness     = Float( 0.03, unit = 'm', input = True, auto_set = False, enter_set = True, table_field = True )
    # age of the concrete at the time of testing
    age           = Int(     28, unit = 'd', input = True, auto_set = False, enter_set = True, table_field = True )
    loading_rate  = Float( 0.50, unit = 'mm/min', input = True, auto_set = False, enter_set = True, table_field = True )

    #--------------------------------------------------------------------------------
    # select the concrete mixture from the concrete database
    #--------------------------------------------------------------------------------

    concrete_mixture_key = Enum( ConcreteMixture.db.keys(), 
                                 input = True, table_field = True )  
    concrete_mixture_ref = Property( Instance( ConcreteMixture ) )
    def _get_concrete_mixture_ref(self):
        return ConcreteMixture.db[ self.concrete_mixture_key ]

    #--------------------------------------------------------------------------------
    # select the textile cross section from the textile database:
    #--------------------------------------------------------------------------------

    fabric_layup_key = Enum( FabricLayup.db.keys(), table_field = True, 
                                      input = True, auto_set = False, enter_set = True  )  
    fabric_layup_ref = Property( Instance( FabricLayup ) )
    def _get_fabric_layup_ref(self ):
        return FabricLayup.db[ self.fabric_layup_key ]

    #--------------------------------------------------------------------------------
    # composite cross section combining the above concrete mixture 
    # and textile cross section:
    #--------------------------------------------------------------------------------

    composite_cross_section_key = Property( Str )
    def _get_composite_cross_section_key(self):
        return self.fabric_layup_key + '_' + self.concrete_mixture_key

    composite_cross_section_ref = Property( Instance( CompositeCrossSection ),
                                            depends_on = 'fabric_layup_key, concrete_mixture_key' )
    def _get_composite_cross_section_ref(self ):
        key = self.composite_cross_section_key
        ref = CompositeCrossSection.db.get( key, None )        
        if ref == None:
            # This material combination was not yet in the material database
            # it gets constructed here but is not persistent yet
            # it will be stored in the database in the __setstate__ method
            # once this test was confirmed for saving 
            return  CompositeCrossSection( fabric_layup_key = self.fabric_layup_key,
                                           concrete_mixture_key = self.concrete_mixture_key )
        else:
            # It was already in the database
            return ref

    #--------------------------------------------------------------------------------
    # calculated material properties 
    #--------------------------------------------------------------------------------
    
    # number of reinforcement layers with orientation in 0-degree direction 
    n_layers_0 = Property(Int, table_field = True)
    def _get_n_layers_0( self ):
        return self.fabric_layup_ref.get_n_layers_0( self.thickness )        

    # number of reinforcement layers with orientation in 90-degree direction 
    n_layers_90 = Property(Int, table_field = True)
    def _get_n_layers_90( self ):
        return self.fabric_layup_ref.get_n_layers_90( self.thickness )        

    # reinforcement ration of the composite 
    rho_c = Property( Float, unit = '-',   depends_on = '+input', table_field = True )
    def _get_rho_c(self):
        return self.composite_cross_section_ref.get_rho_c( self.thickness, self.width )        
    
    # reinforcement ration of the composite (reference value of a periodic cell) 
    rho_cc = Property( Float, unit = '-',   depends_on = '+input', table_field = True )
    def _get_rho_cc(self):
        return self.composite_cross_section_ref.rho_cc        
    
    # E-modulus of the composite at the time of testing 
    E_c = Property( Float, unit = 'MPa', depends_on = '+input', table_field = True  )
    def _get_E_c(self):
        return self.composite_cross_section_ref.get_E_c_time( self.age, self.thickness, self.width )
    
    # E-modulus of the composite after 28 days
    E_c28 = Property( Float, unit = 'MPa', depends_on = '+input', table_field = True  )
    def _get_E_c28(self):
        return self.composite_cross_section_ref.get_E_c_time( 28, self.thickness, self.width )
    
    # E-modulus of the composite at the time of testing 
    E_cc = Property( Float, unit = 'MPa', depends_on = '+input', table_field = True  )
    def _get_E_cc(self):
        return self.composite_cross_section_ref.E_cc

    # E-modulus of the concrete at the time of testing 
    E_m = Property( Float, unit = 'MPa', depends_on = 'age,concrete_mixture_ref' )
    def _get_E_m(self):
        return self.concrete_mixture_ref.get_E_m_time( self.age )
    
    # E-modulus of the concrete  after 28 days 
    E_m28 = Property( Float, unit = 'MPa', depends_on = 'age,concrete_mixture_ref' )
    def _get_E_m28(self):
        return self.concrete_mixture_ref.E_m28
    
    #--------------------------------------------------------------------------------
    # recover the links to the database records
    #--------------------------------------------------------------------------------
    
    def __setstate__(self, dict):
        '''Overload the setstate to recover the links to the database records.
        '''
#        # delete old keys from dict in order to read old pickle files
#        for key_old in ['E_c', 'n_roving_0', 'n_roving_90']:
#            if dict.has_key( key_old ):
#                del dict[ key_old ]
        super( ExPlateTest, self ).__setstate__( dict )

    def __getstate__( self ):
        '''Overload the getstate to recover the links to the database records.
        '''
        dict = super( ExPlateTest, self ).__getstate__()
        key = self.composite_cross_section_key
        ref = CompositeCrossSection.db.get( key, None )
        if ref == None:
            # not in the database - save it now
            CompositeCrossSection.db[ key ] = self.composite_cross_section_ref
        return dict

    #--------------------------------------------------------------------------------
    # define processing
    #--------------------------------------------------------------------------------
     
    @on_trait_change('+input')
    def process_source_data( self ):
        super( ExPlateTest, self ).process_source_data()
        # if center displacement gauge ('WA_M') is missing the measured 
        # displacement of the cylinder ('Weg') is used instead.
        # A minor mistake is made depending on how much time passes
        # before the cylinder has contact with the plate.
        # Verify this in the plot. 
        if 'WA_M' not in self.factor_list:
            print '*** NOTE: Displacement gauge at center ("WA_M") missing. Cylinder displacement ("Weg") is used instead! ***'
#            print 'factor_list', self.factor_list
            self.WA_M = self.Weg

    #--------------------------------------------------------------------------------
    # plot templates
    #--------------------------------------------------------------------------------

    def _plot_templates_default(self):
        return {'force / deflection (center)'          : '_plot_force_deflection_center',
                'smoothed force / deflection (center)' : '_plot_smoothed_force_deflection_center',               
                 }

    def _plot_force_deflection_center(self, axes ):

        xkey = 'deflection [mm]'
        ykey = 'force [kN]'
        xdata = - self.WA_M 
        ydata = - self.Kraft
    
        axes.set_xlabel( '%s' % (xkey,) )
        axes.set_ylabel( '%s' % (ykey,) )
                    
        axes.plot( xdata, ydata
                       # color = c, linewidth = w, linestyle = s 
                       )

    n_fit_window_fraction = Float( 0.1  ) 
            
    def _plot_smoothed_force_deflection_center(self, axes ):
        
        # get the index of the maximum stress
        max_force_idx = argmax( - self.Kraft )
        # get only the ascending branch of the response curve
        f_asc     = - self.Kraft[:max_force_idx + 1]
        w_asc     = - self.WA_M[:max_force_idx + 1]
        
        f_max = f_asc[-1]
        w_max = w_asc[-1]

        n_points = int( self.n_fit_window_fraction * len( w_asc ) )
        f_smooth = smooth( f_asc, n_points, 'flat' )
        w_smooth = smooth( w_asc, n_points, 'flat' )
    
        axes.plot( w_smooth, f_smooth, color = 'blue', linewidth = 2 )

        secant_stiffness_w10 = (f_smooth[10] - f_smooth[0]) / (w_smooth[10] - w_smooth[0] )
        w0_lin = array( [0.0, w_smooth[10] ], dtype = 'float_' )
        f0_lin = array( [0.0, w_smooth[10] * secant_stiffness_w10 ], dtype = 'float_' )

        #axes.plot( w0_lin, f0_lin, color = 'black' )

    #--------------------------------------------------------------------------------
    # view
    #--------------------------------------------------------------------------------

    traits_view = View( VGroup(
                         Group(
                              Item('thickness',   editor=float_editor_3f),
                              Item('edge_length', editor=float_editor_3f),
                              label = 'geometry'
                              ),
                         Group(
                              Item('loading_rate' ),
                              label = 'loading'
                              ),
                         Group(
                              Item('concrete_mixture_key', label = 'concrete mixture'),
                              Item('age'),
                              Item('E_m',   style = 'readonly',  editor=float_editor_0f),
                              Item('E_m28', style = 'readonly',  editor=float_editor_0f),
                              Item('fabric_layup_key', label = 'textile cross section'),
                              VGroup(Item('n_layers_0', style = 'readonly'),
                                     Item('n_layers_90',style = 'readonly')),
                              label = 'material'
                              ),
                         Group(
                              HGroup(Spring(),
                                     Item('composite_cross_section_key', show_label = False, emphasized = True, style = 'readonly'),
                                     Spring()),
                              HGroup(
                              VGroup(
                                     Item('rho_c', style = 'readonly', show_label = True ),
                                     Item('rho_c', style = 'readonly', show_label = True ),
                                     Item('rho_cc',style = 'readonly', show_label = True, editor=float_editor_4f )),
                              VGroup(Item('E_c',   style = 'readonly', show_label = True ),
                                     Item('E_c28', style = 'readonly', show_label = True ),
                                     Item('E_cc',  style = 'readonly', show_label = True, editor=float_editor_0f ))),
                              label = 'composite cross section'
                                  ),                              
                              ),
                        scrollable = True,
                        resizable = True,
                        height = 0.8,
                        width = 0.5,
                        )

if __name__ == '__main__':
    et = ExPlateTest()
    print 'loading_rate', et.loading_rate 
    et.configure_traits()
    