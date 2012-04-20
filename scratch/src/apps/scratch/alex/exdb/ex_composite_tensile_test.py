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
# Created on Feb 15, 2010 by: rch, ascholzen

# @todo - construct the class for fabric layout returning calculating the 
#         cs-area of the reinforcement.
#       - instead of processed array - construct the array traits accessible
#         with the name of the measured channels
#       - reread the pickle file without processing the data (take care to reestablish
#         the link from the ex_type to the ex_run
#       - define the exdb_browser showing the inputs and outputs in a survey
#       - define the ExTreatment class with cummulative evaluation of the response values.
#       
#

from enthought.traits.api import \
    HasTraits, Directory, List, Int, Float, Any, \
    on_trait_change, File, Constant, Instance, Trait, \
    Array, Str, Property, cached_property, WeakRef, \
    Dict, Button, Bool, Enum, Event, implements, \
    DelegatesTo

from enthought.util.home_directory import \
    get_home_directory

from enthought.traits.ui.api import \
    View, Item, DirectoryEditor, TabularEditor, HSplit, VGroup, \
    TableEditor, EnumEditor, Handler, FileEditor, VSplit, Group, \
    InstanceEditor, HGroup, Spring

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

from numpy import \
    array, fabs, where, copy, ones, linspace, ones_like, hstack

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
from i_ex_type import IExType

from mathkit.array.smoothing import smooth

from promod.matdb.trc.fabric_layup \
    import FabricLayup
    
from promod.matdb.trc.fabric_layout \
    import FabricLayout 
  
from promod.matdb.trc.concrete_mixture \
    import ConcreteMixture
    
from promod.matdb.trc.composite_cross_section \
    import CompositeCrossSection
  
    
class ExCompositeTensileTest( ExType ):
    '''Read the data from the directory
    '''
    
    implements( IExType )

    #--------------------------------------------------------------------------------
    # specify inputs:
    #--------------------------------------------------------------------------------

    width         = Float( 0.100, unit = 'm', table_field = True,
                           input = True, auto_set = False, enter_set = True )
    thickness     = Float( 0.030, unit = 'm', table_field = True,
                           input = True, auto_set = False, enter_set = True )
    gauge_length  = Float( 0.450, unit = 'm', table_field = True,
                           input = True, auto_set = False, enter_set = True )
    
    # age of the concrete at the time of testing
    age           = Int(    28, unit = 'd',  table_field = True,
                            input = True, auto_set = False, enter_set = True )
    loading_rate  = Float( 1.0, unit =  'mm/min', table_field = True,
                           input = True, auto_set = False, enter_set = True )

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

    fabric_layup_key = Enum( FabricLayup.db.keys(), 
                                      input = True, table_field = True )  
    fabric_layup_ref = Property( Instance( FabricLayup ) )
    def _get_fabric_layup_ref(self ):
        return FabricLayup.db[ self.fabric_layup_key ]

    # textile fabric layout used to build the textile cross section
    #
    fabric_layout = Property( Instance(FabricLayout), 
                                      depends_on = 'fabric_layup_key' )
    def _get_fabric_layout( self ):
        return self.fabric_layup_ref.fabric_layout
    
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

    #--------------------------------------------------------------------------
    # Get properties for the specified width and thickness
    #--------------------------------------------------------------------------

    # Discrete number of rovings used in the tensile specimens in 
    # longitudinal direction for ONE LAYER of textile fabric of size 'width'.
    # The rovings are distiguished in 0- and  90-direction, respectively. 
    # The value is calculated based on spacing of the rovings specified in
    # 'textil_fabric_layout' and the defined 'width' rounded to the next smaller
    # integer. 
    #
    
    n_rovings_0  = Property( Int, unit = '-', table_field = True,
                             input = True, auto_set = False, enter_set = True ) 
    def _get_n_rovings_0( self ):
        return self.fabric_layout.get_n_rovings_0( self.width )

    n_rovings_90  = Property( Int, unit = '-', table_field = True,
                              input = True, auto_set = False, enter_set = True ) 
    def _get_n_rovings_90( self ):
        return self.fabric_layout.get_n_rovings_90( self.width )        
    
    # number of reinforcement layers with orientation in 0-degree direction 
    n_layers_0 = Property(Int, unit = '-', table_field = True)
    def _get_n_layers_0( self ):
        return self.fabric_layup_ref.get_n_layers_0( self.thickness )        

    # number of reinforcement layers with orientation in 90-degree direction 
    n_layers_90 = Property(Int, unit = '-', table_field = True)
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
    
    # cross-sectional-area of the composite 
    A_c = Property( Float, unit = 'm^2' )
    def _get_A_c(self):
        return self.width * self.thickness 
        
    # total cross-sectional-area of the textile reinforcement 
    A_tex = Property( Float, unit = 'mm^2' )
    def _get_A_tex(self):
        return self.composite_cross_section_ref.get_A_tex( self.thickness, self.width  )
    
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
        super( ExCompositeTensileTest, self ).__setstate__( dict )

    def __getstate__( self ):
        '''Overload the getstate to recover the links to the database records.
        '''
        dict = super( ExCompositeTensileTest, self ).__getstate__()
        key = self.composite_cross_section_key
        ref = CompositeCrossSection.db.get( key, None )
        if ref == None:
            # not in the database - save it now
            CompositeCrossSection.db[ key ] = self.composite_cross_section_ref
        return dict

    #--------------------------------------------------------------------------------
    # Output characteristics
    #--------------------------------------------------------------------------------

    sig_c_max   = Float( 0.0, output = True, table_field = True, unit = 'MPa' )
    sig_tex_max = Float( 0.0, output = True, table_field = True, unit = 'MPa' )
    eps_c_max   = Float( 0.0, output = True, table_field = True, unit = '-' )

    #--------------------------------------------------------------------------------
    # define processing
    #--------------------------------------------------------------------------------

    @on_trait_change('+input')
    def process_source_data( self ):
        super( ExCompositeTensileTest, self ).process_source_data()
        # 
        self.W10_re -= self.W10_re[0]
        self.W10_re *= -1
        
        self.W10_li -= self.W10_li[0]
        self.W10_li *= -1

        self.W10_vo -= self.W10_vo[0]
        self.W10_vo *= -1

        self._average_strain()
        self._sig_c()
        self._sig_tex()
        self._get_ascending_branch()
        self._get_polyfit()

    eps = Array( 'float_', output = True )
    def _average_strain( self ):
        # measured strains 
        eps_li = self.W10_li / ( self.gauge_length * 1000. )  #[mm/mm]
        eps_re = self.W10_re / ( self.gauge_length * 1000. )
        eps_vo = self.W10_vo / ( self.gauge_length * 1000. )
        
        # average strains 
        eps_m = ( (eps_li + eps_re)/2. + eps_vo)/2.
        min_eps = min( eps_m[:10] )
        self.eps = eps_m - min_eps        

    sig_c = Array( 'float_', output = True )
    def _sig_c( self ):
        # measured force: 
        force = self.Kraft # [kN]
        # cross sectional area of the concrete [m^2]: 
        A_c = self.A_c 
        # calculated stress: 
        sig_c = ( force / 1000. ) / A_c  # [MPa]
        self.sig_c = sig_c        

    sig_tex = Array( 'float_', output = True )
    def _sig_tex( self ):
        # measured force: 
        force = self.Kraft # [kN]
        # cross sectional area of the reinforcement: 
        A_tex = self.A_tex 
        # calculated stresses:
        sig_tex = ( force * 1000. ) / self.A_tex  # [MPa]
        self.sig_tex = sig_tex
                
    eps_asc     = Array( 'float_' )
    sig_c_asc   = Array( 'float_', output = True )
    sig_tex_asc = Array( 'float_', output = True )
    def _get_ascending_branch( self ):
        # get the index of the maximum stress
        max_stress_idx = argmax( self.sig_c )
        # get only the ascending branch of the response curve
        self.eps_asc       = self.eps[:max_stress_idx + 1]
        self.sig_c_asc     = self.sig_c[:max_stress_idx + 1]
        self.sig_tex_asc   = self.sig_tex[:max_stress_idx + 1]
        
        self.eps_c_max = self.eps_asc[-1]
        self.sig_tex_max = self.sig_tex_asc[-1]
        self.sig_c_max = self.sig_c_asc[-1]

    eps_fit     = Array( 'float_', output = True )
    sig_c_fit   = Array( 'float_', output = True )
    sig_tex_fit = Array( 'float_', output = True )
    n_poly = Int(12)
    n_fit_points = Int( 100 )
    n_fit_window_fraction = Float( 0.1  ) 
    def _get_polyfit( self ):
        # get the fit with n-th-order polynomial
        
        n_points = int( self.n_fit_window_fraction * len( self.eps ) )

        self.eps_fit = smooth( self.eps_asc,   n_points, 'flat' )
        sig_c_fit    = smooth( self.sig_c_asc, n_points, 'flat' )

        sig_lin = self.E_c * self.eps_fit
        cut_sig = where( sig_c_fit > sig_lin )
        sig_c_fit[ cut_sig ] = sig_lin[ cut_sig ]
        self.sig_c_fit = sig_c_fit
        
        self.sig_tex_fit = smooth( self.sig_tex_asc, n_points, 'flat' )

    #--------------------------------------------------------------------------------
    # plot templates
    #--------------------------------------------------------------------------------

    plot_templates = {'concrete stress / strain'           : '_plot_c_stress_strain',
                      'smoothed concrete stress / strain'  : '_plot_c_fit_stress_strain',               
                      'textile stress / strain'            : '_plot_tex_stress_strain',               
                      'smoothed textile stress / strain'   : '_plot_tex_fit_stress_strain',               
                      'force / gauge displacement'         : '_plot_force_displacement' }
    
    def _plot_c_stress_strain(self, axes ):

        xkey = 'eps [-]'
        ykey = 'sig_c [MPa]'
        xdata = self.eps 
        ydata = self.sig_c
    
        axes.set_xlabel( '%s' % (xkey,) )
        axes.set_ylabel( '%s' % (ykey,) )
                    
        axes.plot( xdata, ydata
                       # color = c, linewidth = w, linestyle = s 
                       )

    def _plot_c_fit_stress_strain(self, axes ):

        #self._plot_c_stress_strain(axes)

        axes.set_xlabel( 'eps_asc [-]' )
        axes.set_ylabel( 'sig_c_fit [MPa]')

        axes.plot( self.eps_asc, self.sig_c_asc, color = 'green' 
                       # color = c, linewidth = w, linestyle = s 
                       )                      
        sig_lin = array([0, self.sig_c_max], dtype = 'float_' )
        eps_lin = array([0, self.sig_c_max / self.E_c ], dtype = 'float_' )
        axes.plot( eps_lin, sig_lin, color = 'red' )
                    
        axes.plot( self.eps_fit, self.sig_c_fit, color = 'blue', linewidth = 2 
                       # color = c, linewidth = w, linestyle = s 
                       )

    def _plot_tex_stress_strain(self, axes ):

        xkey = 'eps [-]'
        ykey = 'sig_tex [MPa]'
        xdata = self.eps
        ydata = self.sig_tex
    
        axes.set_xlabel( '%s' % (xkey,) )
        axes.set_ylabel( '%s' % (ykey,) )
                    
        axes.plot( xdata, ydata
                       # color = c, linewidth = w, linestyle = s 
                       )

    def _plot_tex_fit_stress_strain(self, axes ):

        axes.set_xlabel( 'eps_asc [-]' )
        axes.set_ylabel( 'sig_tex_fit [MPa]')

        axes.plot( self.eps_asc, self.sig_tex_asc
                       # color = c, linewidth = w, linestyle = s 
                       )
                    
        axes.plot( self.eps_fit, self.sig_tex_fit
                       # color = c, linewidth = w, linestyle = s 
                       )
        
        eps_lin = array( [0, self.eps_fit[-1]], dtype = 'float_' )
        sig_lin = self.fabric_layout.E_tex_0 * eps_lin
        # plot the textile secant stiffness at fracture state 
        axes.plot( eps_lin, sig_lin
                       # color = c, linewidth = w, linestyle = s 
                       )

    def _plot_force_displacement(self, axes ):        

        # 
        axes.plot( self.W10_re, self.Kraft )
        axes.plot( self.W10_li, self.Kraft )
        axes.plot( self.W10_vo, self.Kraft )
        axes.set_xlabel( '%s' % ('Weg [mm]',) )
        axes.set_ylabel( '%s' % ('Kraft [kN]',) )
            
    traits_view = View( VSplit(
                         Group(
                               Group(
                                  Item('width',        editor=float_editor_3f),
                                  Item('thickness',    editor=float_editor_3f),
                                  Item('gauge_length', editor=float_editor_3f),
                                  label = 'geometry'
                                  ),
                               Group(
                                  Item('loading_rate'),
                                  label = 'loading'
                                  ),
                               Group(
                                  Item('concrete_mixture_key', label = 'concrete mixture'),
                                  Item('age'),
                                  Item('E_m',   style = 'readonly', editor=float_editor_0f ),
                                  Item('E_m28', style = 'readonly', editor=float_editor_0f),
                                  Item('fabric_layup_key', label = 'textile cross section'),
                                  HGroup(VGroup(Item('n_layers_0',   style = 'readonly'),
                                                Item('n_layers_90',  style = 'readonly')),
                                         VGroup(Item('n_rovings_0',  style = 'readonly'),
                                                Item('n_rovings_90', style = 'readonly'))),
                                  label = 'material'
                                  ),
                               Group(
                                  HGroup(Spring(),
                                     Item('composite_cross_section_key', show_label = False, emphasized = True, style = 'readonly'),
                                     Spring()),     
                                  HGroup(
                                  VGroup(Item('rho_c',  style = 'readonly', show_label = True, editor=float_editor_4f  ),
                                         Item('rho_cc', style = 'readonly', show_label = True, editor=float_editor_4f )),
                                  VGroup(Item('E_c',    style = 'readonly', show_label = True, editor=float_editor_0f ),
                                         Item('E_c28',  style = 'readonly', show_label = True, editor=float_editor_0f ),
                                         Item('E_cc',   style = 'readonly', show_label = True, editor=float_editor_0f ))),
                                  label = 'composite cross section'
                                  ),                              
                               label = 'input variables',
                               id = 'promod.exdb.ex_composite_tensile_test.vgroup.inputs',
                               dock = 'tab',
                               scrollable = True,
                               ),
                         Group(
                               Item('sig_c_max',   
                                   emphasized = True, style = 'readonly', editor=float_editor_2f ),
                               Item('sig_tex_max', 
                                   emphasized = True, style = 'readonly', editor=float_editor_2f ),
                               Item('eps_c_max',      
                                   emphasized = True, style = 'readonly', editor=float_editor_4f  ),
                               Item('E_c',                
                                   emphasized = True, style = 'readonly', editor=float_editor_0f  ),
                               label = 'output characteristics',
                               id = 'promod.exdb.ex_composite_tensile_test.vgroup.outputs',
                               dock = 'tab',
                               scrollable = True,
                         ),
                         scrollable = True,
                         id = 'promod.exdb.ex_composite_tensile_test.vgroup',
                         dock = 'tab',
                         ),
                         id = 'promod.exdb.ex_composite_tensile_test',
                         dock = 'tab',
                         scrollable = True,
                         resizable = True,
                         height = 0.8,
                         width = 0.5,
                         )
    
if __name__ == '__main__':
    et = ExCompositeTensileTest()
    print 'E_c', et.E_c 
    print 'fabric_layout', et.fabric_layout
    et.configure_traits()
    