'''
Created on Jun 23, 2010

@author: alexander
'''

from enthought.traits.api import \
    HasTraits, Directory, List, Int, Float, Any, Enum, \
    on_trait_change, File, Constant, Instance, Trait, \
    Array, Str, Property, cached_property, WeakRef, \
    Dict, Button, Color, Bool, Trait

from enthought.util.home_directory import \
    get_home_directory

from enthought.traits.ui.api import \
    View, Item, DirectoryEditor, TabularEditor, HSplit, Tabbed, VGroup, \
    TableEditor, Group, ListEditor, VSplit, HSplit, VGroup, HGroup, Spring, \
    Include

from enthought.mayavi import \
    mlab

#from enthought.mayavi.mlab import \
#    colorbar, show, points3d
#
#from enthought.mayavi.api import \
#    Engine


from enthought.traits.ui.table_column import \
    ObjectColumn

from enthought.traits.ui.menu import \
    OKButton, CancelButton

from enthought.traits.ui.tabular_adapter \
    import TabularAdapter

from numpy import array, loadtxt, arange, sqrt, zeros, arctan, sin, cos, ones_like, \
                vstack, savetxt, hstack, argsort, fromstring, zeros_like, \
                copy, c_, newaxis, argmax, where, argsort, sqrt, frompyfunc, \
                max as ndmax, ones

from math import pi
from string import split
import os

from scipy.io import read_array

from get_sig_fl_iterative_rch import SigFlCalib

class LSArrayAdapter ( TabularAdapter ):

    columns = Property
    def _get_columns( self ):
#        print 'GETTING COLUMNS', self.object.columns, self.object, self.object.__class__
        columns = self.object.columns
        return [ ( name, idx ) for idx, name in enumerate( columns ) ]

    font = 'Courier 10'
    alignment = 'right'
    format = '%5.2f'#'%g'
    even_bg_color = Color( 0xE0E0FF )
    width = Float( 80 )

    #@todo: format columns using 'column_id'
#    adapter_column_map = Property(depends_on = 'adapters,columns')


class LS( HasTraits ):
    '''Limit state class
    '''

    # backward link to the info shell to access the
    # input data when calculating 
    # the limit-state-specific values
    #
    ls_table = WeakRef

    #-------------------------------
    # ls columns
    #-------------------------------
    # defined in the subclasses
    #
    ls_columns = List
    show_ls_columns = Bool( True )

    #-------------------------------
    # geo columns form info shell
    #-------------------------------

    geo_columns = List( [ 'elem_no', 'X', 'Y', 'Z', 'thickness' ] )
    show_geo_columns = Bool( True )

    elem_no = Property( Array )
    def _get_elem_no( self ):
        return self.ls_table.elem_no

    X = Property( Array )
    def _get_X( self ):
        return self.ls_table.X

    Y = Property( Array )
    def _get_Y( self ):
        return self.ls_table.Y

    Z = Property( Array )
    def _get_Z( self ):
        return self.ls_table.Z

    thickness = Property( Array )
    def _get_thickness( self ):
        return self.ls_table.thickness

    #-------------------------------
    # state columns form info shell
    #-------------------------------

    state_columns = List( [
#                           'mx', 'my', 'mxy', 'nx', 'ny', 'nxy',
#                           'sigx_lo', 'sigy_lo', 'sigxy_lo', 
                           'sig1_lo', 'sig1_up_sig_lo', 'alpha_sig_lo',
                           'm_sig_lo', 'n_sig_lo',
#                           'sigx_up', 'sigy_up', 'sigxy_up', 
                           'sig1_up', 'sig1_lo_sig_up', 'alpha_sig_up',
                           'm_sig_up', 'n_sig_up',
                            ] )

    show_state_columns = Bool( True )

    mx = Property( Array )
    def _get_mx( self ):
        return self.ls_table.mx

    my = Property( Array )
    def _get_my( self ):
        return self.ls_table.my

    mxy = Property( Array )
    def _get_mxy( self ):
        return self.ls_table.mxy

    nx = Property( Array )
    def _get_nx( self ):
        return self.ls_table.nx

    ny = Property( Array )
    def _get_ny( self ):
        return self.ls_table.ny

    nxy = Property( Array )
    def _get_nxy( self ):
        return self.ls_table.nxy


    # evaluate principal stresses
    # upper face:
    #
    sigx_up = Property( Array )
    def _get_sigx_up( self ):
        return self.ls_table.sigx_up

    sigy_up = Property( Array )
    def _get_sigy_up( self ):
        return self.ls_table.sigy_up

    sigxy_up = Property( Array )
    def _get_sigxy_up( self ):
        return self.ls_table.sigxy_up

    sig1_up = Property( Array )
    def _get_sig1_up( self ):
        return self.ls_table.sig1_up

    sig2_up = Property( Array )
    def _get_sig2_up( self ):
        return self.ls_table.sig2_up

    alpha_sig_up = Property( Array )
    def _get_alpha_sig_up( self ):
        return self.ls_table.alpha_sig_up

    # lower face:
    #
    sigx_lo = Property( Float )
    def _get_sigx_lo( self ):
        return self.ls_table.sigx_lo

    sigy_lo = Property( Float )
    def _get_sigy_lo( self ):
        return self.ls_table.sigy_lo

    sigxy_lo = Property( Float )
    def _get_sigxy_lo( self ):
        return self.ls_table.sigxy_lo

    sig1_lo = Property( Float )
    def _get_sig1_lo( self ):
        return self.ls_table.sig1_lo

    sig2_lo = Property( Float )
    def _get_sig2_lo( self ):
        return self.ls_table.sig2_lo

    alpha_sig_lo = Property( Float )
    def _get_alpha_sig_lo( self ):
        return self.ls_table.alpha_sig_lo

    #-------------------------------
    # ls table
    #-------------------------------

    # all columns associated with the limit state including the corresponding
    # stress resultants
    #
    columns = Property( List, depends_on = 'show_geo_columns, show_state_columns, show_ls_columns' )
    @cached_property
    def _get_columns( self ):
        columns = []

        if self.show_geo_columns:
            columns += self.geo_columns

        if self.show_state_columns:
            columns += self.state_columns

        if self.show_ls_columns:
            columns += self.ls_columns

        return columns

    # select column used for sorting the data in selected sorting order 
    #
    sort_column = Enum( values = 'columns' )
    def _sort_column_default( self ):
        return self.columns[-1]

    sort_order = Enum( 'descending', 'ascending', 'unsorted' )

    #-------------------------------------------------------
    # get the maximum value of the selected variable 
    # 'max_in_column' of the current sheet (only one sheet)
    #-------------------------------------------------------

    # get the maximum value of the chosen column
    #
    max_in_column = Enum( values = 'columns' )
    def _max_in_column_default( self ):
        return self.columns[-1]

    max_value = Property( depends_on = 'max_in_column' )
    def _get_max_value( self ):
        col = getattr( self, self.max_in_column )[:, 0]
        return max( col )

    #-------------------------------------------------------
    # get the maximum value and the corresponding case of 
    # the selected variable 'max_in_column' in all (!) sheets
    #-------------------------------------------------------

    max_value_all = Property( depends_on = 'max_in_column' )
    def _get_max_value_all( self ):
        return self.ls_table.max_value_and_case[ self.max_in_column ]['max_value']

    max_case = Property( depends_on = 'max_in_column' )
    def _get_max_case( self ):
        return self.ls_table.max_value_and_case[ self.max_in_column ]['max_case']

    #-------------------------------------------------------
    # get ls_table for View
    #-------------------------------------------------------

    # stack columns together for table used by TabularEditor
    #
    ls_array = Property( Array, depends_on = 'sort_column, sort_order, \
                                              show_geo_columns, \
                                              show_state_columns, \
                                              show_ls_columns' )

    @cached_property
    def _get_ls_array( self ):

        arr_list = [ getattr( self, col ) for col in self.columns ]

        # get the array currently selected by the sort_column enumeration
        #
        sort_arr = getattr( self, self.sort_column )[:, 0]
        sort_idx = argsort( sort_arr )
        ls_array = hstack( arr_list )

        if self.sort_order == 'descending':
            return ls_array[ sort_idx[::-1] ]
        if self.sort_order == 'ascending':
            return ls_array[ sort_idx ]
        if self.sort_order == 'unsorted':
            return ls_array

    #---------------------------------
    # plot outputs in mlab-window 
    #---------------------------------

    plot_column = Enum( values = 'columns' )
    plot = Button
    def _plot_fired( self ):
        X = self.ls_table.X[:, 0]
        Y = self.ls_table.Y[:, 0]
        Z = self.ls_table.Z[:, 0]
        plot_col = getattr( self, self.plot_column )[:, 0]

        if self.plot_column == 'n_tex':
            plot_col = where( plot_col < 0, 0, plot_col )

        mlab.figure( figure = "SFB532Demo",
                     bgcolor = ( 1.0, 1.0, 1.0 ),
                     fgcolor = ( 0.0, 0.0, 0.0 ) )

        mlab.points3d( X, Y, ( -1.0 ) * Z, plot_col,
#                       colormap = "gist_rainbow",
#                       colormap = "Reds",
                       colormap = "YlOrBr",
                       mode = "cube",
                       scale_factor = 0.15 )

        mlab.scalarbar( title = self.plot_column, orientation = 'vertical' )

        mlab.show



    # name of the trait that is used to assess the evaluated design
    #
    assess_name = Str( '' )

    #-------------------------------
    # ls group
    #-------------------------------

    # @todo: the dynamic selection of the columns to be displayed 
    # does not work in connection with the LSArrayAdapter 
    ls_group = VGroup( 
                        HGroup( #Item( 'assess_name' ),
                                Item( 'max_in_column' ),
                                Item( 'max_value', style = 'readonly', format_str = '%6.2f' ),
                              ),
                        HGroup( Item( 'sort_column' ),
                                Item( 'sort_order' ),
                                Item( 'show_geo_columns', label = 'show geo' ),
                                Item( 'show_state_columns', label = 'show state' ),
                                Item( 'show_ls_columns', label = 'show ls' ),
                                Item( 'plot_column' ),
                                Item( 'plot' ),
                              ),
                     )


class SLS( LS ):
    '''Serviceability limit state
    '''

    # ------------------------------------------------------------
    # SLS: material parameters (Inputs)
    # ------------------------------------------------------------

    # tensile strength [MPa]
    f_ctk = Float( 4.0, input = True )

    # flexural tensile strength [MPa]
    f_m = Float( 5.0, input = True )

    # ------------------------------------------------------------
    # SLS - derived params:
    # ------------------------------------------------------------

    # area
    #
    A = Property( Float )
    def _get_A( self ):
        return self.ls_table.thickness * 1.

    # moment of inertia
    #
    W = Property( Float )
    def _get_W( self ):
        return 1. * self.ls_table.thickness ** 2 / 6.

    # ------------------------------------------------------------
    # SLS: outputs
    # ------------------------------------------------------------

    ls_columns = List( ['sig_n', 'sig_m', 'eta_n', 'eta_m', 'eta_tot', ] )

    ls_values = Property( depends_on = '+input' )
    @cached_property
    def _get_ls_values( self ):
        '''get the outputs for SLS
        '''
        n = self.n
        m = self.m

        A = self.A
        W = self.W
        f_ctk = self.f_ctk
        f_m = self.f_m

        sig_n = n / A / 1000.
        sig_m = abs( m / W ) / 1000.
        eta_n = sig_n / f_ctk
        eta_m = sig_m / f_m
        eta_tot = eta_n + eta_m

        return { 'sig_n':sig_n, 'sig_m':sig_m,
                 'eta_n':eta_n, 'eta_m':eta_m,
                 'eta_tot':eta_tot }

    sig_n = Property
    def _get_sig_n( self ):
        return self.ls_values['sig_n']

    sig_m = Property
    def _get_sig_m( self ):
        return self.ls_values['sig_m']

    eta_n = Property
    def _get_eta_n( self ):
        return self.ls_values['eta_n']

    eta_m = Property
    def _get_eta_m( self ):
        return self.ls_values['eta_m']

    eta_tot = Property
    def _get_eta_tot( self ):
        return self.ls_values['eta_tot']

    assess_name = 'max_eta_tot'

#    assess_name = 'max_sig1_up'
#
#    max_sig1_up = Property( depends_on = '+input' )
#    @cached_property
#    def _get_max_sig1_up( self ):
#        return ndmax( self.sig1_up )

    # @todo: make it possible to select the assess value:
    #
#    assess_name = Enum( values = 'columns' )
#    def _assess_name_default( self ):
#        return self.columns[-1]

    max_eta_tot = Property( depends_on = '+input' )
    @cached_property
    def _get_max_eta_tot( self ):
        return ndmax( self.eta_tot )


    #-------------------------------
    # ls view
    #-------------------------------

    # @todo: the dynamic selection of the columns to be displayed 
    # does not work in connection with the LSArrayAdapter 
    traits_view = View( VGroup( 
                            HGroup( Item( name = 'f_ctk', label = 'Tensile strength concrete [MPa]: f_ctk ' ),
                                    Item( name = 'f_m', label = 'Flexural tensile trength concrete [MPa]: f_m ' )
                                   ),
                            VGroup( 
                                Include( 'ls_group' ),

                                # @todo: currently LSArrayAdapter must be called both 
                                #        in SLS and ULS separately to configure columns 
                                #        arrangement individually
                                #
                                Item( 'ls_array', show_label = False,
                                      editor = TabularEditor( adapter = LSArrayAdapter() ) )
                                  ),
                              ),
                      resizable = True,
                      scrollable = True,
                      height = 1000,
                      width = 1100
                      )

class ULS( LS ):
    '''Ultimate limit state
    '''

    #--------------------------------------------------------
    # ULS: material parameters (Inputs)
    #--------------------------------------------------------

    # tensile composite strength in 0-direction [MPa]
    # derived applying EN-DIN 1990
    # F_t,exp,m = 103.4 kN # mean value tensile test
    # F_Rk,0 = 86.6 kN # characteristic value
    # F_Rd,0 = 86.6 kN / 1.5 = 57.7 kN # design value
    # sig_Rd,0,Zug = 57.7 / 14cm / 6cm = 6.87 MPa
    #
    sig_comp_0_Rd = Float( 6.87 )

    # tensile composite strength in 90-direction [MPa]
    # use value of the 0-direction as conservative simplification:
    #
    sig_comp_90_Rd = Float( 6.87 )

    # tensile strength of the textile reinforcement [kN/m]
    # --> f_Rtex,0 = 6.87 MPa / 10. * 100cm * 6 cm / 12 layers = 34.3 kN/m
    #
    f_Rtex_0 = f_Rtex_90 = 34.3
    print 'NOTE: f_Rtex_0 = f_Rtex_90 = set to %g kN/m !' % ( f_Rtex_0 )

    # sig_Rd,0,Biegung = 7.95 MPa
    k_fl = 7.95 / 6.87
    print 'NOTE: k_fl = set to %g [-] !' % ( k_fl )

    # ------------------------------------------------------------
    # ULS - derived params:
    # ------------------------------------------------------------

    # Parameters for the cracked state (GdT):
    # assumptions!

    # (resultierende statische Nutzhoehe) 
    #
    d = Property( Float )
    def _get_d( self ):
        return 0.75 * self.ls_table.thickness

    # distance from the center of gravity to the resulting reinforcement layer
    #
    zs = Property( Float )
    def _get_zs( self ):
        return self.d - self.ls_table.thickness / 2.

    # inner cantilever
    #
    z = Property( Float )
    def _get_z( self ):
        return 0.9 * self.d

    # ------------------------------------------------------------
    # ULS: outputs
    # ------------------------------------------------------------

    # sig1_lo -direction:
    #
    sig1_up_sig_lo = Property( Array )
    def _get_sig1_up_sig_lo( self ):
        return self.ls_table.sig1_up_sig_lo

    m_sig_lo = Property( Array )
    def _get_m_sig_lo( self ):
        return self.ls_table.m_sig_lo

    n_sig_lo = Property( Array )
    def _get_n_sig_lo( self ):
        return self.ls_table.n_sig_lo

    # sig1_up -direction:
    #
    sig1_lo_sig_up = Property( Array )
    def _get_sig1_lo_sig_up( self ):
        return self.ls_table.sig1_lo_sig_up

    m_sig_up = Property( Array )
    def _get_m_sig_up( self ):
        return self.ls_table.m_sig_up

    n_sig_up = Property( Array )
    def _get_n_sig_up( self ):
        return self.ls_table.n_sig_up

    # choose evaluation mode to calculate the number of reinf-layers 'n_tex':
    #
#    eval_mode = 'massivbau'
    eval_mode = 'princ_sig_level_1'
#    eval_mode = 'princ_sig_level_2'


    ls_values = Property( depends_on = '+input' )
    @cached_property
    def _get_ls_values( self ):
        '''get the outputs for ULS
        '''

        #---------------------------------------------------------
        # conditions for case distinction
        # (-- pure tension -- bending -- compression --) 
        # NOTE: based in all cases of 'eval_mode' on the direction of 
        #       the maximum principle (tensile) stress
        #---------------------------------------------------------

        #----------------------------------------        
        # upper side:
        #----------------------------------------        

        # NOTE: the case zero stresses at both sides would be included more then once
        # irrelevant for real situations
        cond_s1u_ge_0 = self.sig1_up >= 0. # tensile stress upper side
        cond_s1u_le_0 = self.sig1_up <= 0. # compressive stress upper side
        cond_sl_ge_0 = self.sig1_lo_sig_up >= 0. # corresponding tensile stress lower side
        cond_sl_le_0 = self.sig1_lo_sig_up <= 0. # corresponding compressive stress lower side
        cond_s1u_gt_sl = abs( self.sig1_up ) > abs( self.sig1_lo_sig_up )
        cond_s1u_lt_sl = abs( self.sig1_up ) < abs( self.sig1_lo_sig_up )
        cond_s1u_eq_0 = self.sig1_up == 0.
        cond_sl_eq_0 = self.sig1_lo_sig_up == 0.

        # consider only cases for bending without zero stresses at both sides:
        #----------------------------------------        
        cond_up_eq_0 = cond_s1u_eq_0 * cond_sl_eq_0
        cond_up_neq_0 = cond_up_eq_0 - True

        # tension upper side (sig1_up >= 0):
        #----------------------------------------        
        #
        # caused by pure tension:
        #
        cond_up_t = cond_s1u_ge_0 * cond_sl_ge_0
        #
        # caused by bending:
        #
        cond_up_tb = cond_s1u_ge_0 * cond_sl_le_0 * cond_up_neq_0
        #
        # caused by bending and a normal tension:
        #
        cond_up_tb_t = cond_up_tb * cond_s1u_gt_sl * cond_up_neq_0
        #
        # caused by bending and a normal compression:
        #
        cond_up_tb_c = cond_up_tb * cond_s1u_lt_sl * cond_up_neq_0

        # compression upper side (sig1_up <= 0):
        #----------------------------------------        
        #
        # caused by pure compression:
        #
        cond_up_c = cond_s1u_le_0 * cond_sl_le_0
        #
        # caused by bending:
        #
        cond_up_cb = cond_s1u_le_0 * cond_sl_ge_0


        #----------------------------------------        
        # lower side:
        #----------------------------------------        

        # NOTE: the case zero stresses at both sides would be included more then once
        # irrelevant for real situations
        cond_s1l_ge_0 = self.sig1_lo >= 0. # tensile stress lower side
        cond_s1l_le_0 = self.sig1_lo <= 0. # compressive stress lower side
        cond_su_ge_0 = self.sig1_up_sig_lo >= 0. # corresponding tensile stress upper side
        cond_su_le_0 = self.sig1_up_sig_lo <= 0. # corresponding compressive stress upper side
        cond_s1l_gt_su = abs( self.sig1_lo ) > abs( self.sig1_up_sig_lo )
        cond_s1l_lt_su = abs( self.sig1_lo ) < abs( self.sig1_up_sig_lo )
        cond_s1l_eq_0 = self.sig1_lo == 0. # zero stresses at lower side
        cond_su_eq_0 = self.sig1_up_sig_lo == 0.# corresponding zero stresses at upper side

        # consider only cases for bending without zero stresses at both sides:
        #----------------------------------------        
        cond_lo_eq_0 = cond_s1l_eq_0 * cond_su_eq_0
        cond_lo_neq_0 = cond_lo_eq_0 - True

        # tension lower side (sig1_lo >= 0):
        #----------------------------------------        
        #
        # caused by pure tension:
        #
        cond_lo_t = cond_s1l_ge_0 * cond_su_ge_0
        #
        # caused by bending:
        #
        cond_lo_tb = cond_s1l_ge_0 * cond_su_le_0 * cond_lo_neq_0
        #
        # caused by bending and a normal tension:
        #
        cond_lo_tb_t = cond_lo_tb * cond_s1l_gt_su * cond_lo_neq_0
        #
        # caused by bending and a normal compression:
        #
        cond_lo_tb_c = cond_lo_tb * cond_s1l_lt_su * cond_lo_neq_0

        # compression lower side (sig1_lo <= 0):
        #----------------------------------------        
        #
        # caused by pure compression:
        #
        cond_lo_c = cond_s1l_le_0 * cond_su_le_0
        #
        # caused by bending:
        #
        cond_lo_cb = cond_s1l_le_0 * cond_su_ge_0


        #----------------------------------------        
        # check if all elements are classified in one of the cases
        # 'bending, compression, tension' for the upper and the lower side
        # sum of all conditions must be equal to n_elems * 2 (for upper and lower side)
        #----------------------------------------        
#        
#        print 'self.sig1_lo.shape[0]', self.sig1_lo.shape[0]
#
#        print 'cond_lo_t', self.sig1_lo[cond_lo_t].shape[0]
#        print 'cond_lo_tb_t', self.sig1_lo[cond_lo_tb_t].shape[0]
#        print 'cond_lo_tb_c', self.sig1_lo[cond_lo_tb_c].shape[0]
#        print 'cond_lo_tb', self.sig1_lo[cond_lo_tb].shape[0]
#        print 'cond_lo_c', self.sig1_lo[cond_lo_c].shape[0]
#        print 'cond_lo_cb', self.sig1_lo[cond_lo_cb].shape[0]
#
#        print 'cond_up_t', self.sig1_lo[cond_up_t].shape[0]
#        print 'cond_up_tb_t', self.sig1_lo[cond_up_tb_t].shape[0]
#        print 'cond_up_tb_c', self.sig1_lo[cond_up_tb_c].shape[0]
#        print 'cond_up_tb', self.sig1_lo[cond_up_tb].shape[0]
#        print 'cond_up_c', self.sig1_lo[cond_up_c].shape[0]
#        print 'cond_up_cb', self.sig1_lo[cond_up_cb].shape[0]
#
#        print 'sum_lo', \
#                      self.sig1_lo[cond_lo_t].shape[0] + \
#                      self.sig1_lo[cond_lo_tb].shape[0] + \
#                      self.sig1_lo[cond_lo_c].shape[0] + \
#                      self.sig1_lo[cond_lo_cb].shape[0]
#
#        print 'sum_up', \
#                      self.sig1_lo[cond_up_t].shape[0] + \
#                      self.sig1_lo[cond_up_tb].shape[0] + \
#                      self.sig1_lo[cond_up_c].shape[0] + \
#                      self.sig1_lo[cond_up_cb].shape[0]

        #---------------------------------------------------------
        # initialize arrays to be filled by case distinction:
        #---------------------------------------------------------
        #
        f_t_sig_up = zeros_like ( self.sig1_up ) # [kN/m]
        f_t_sig_lo = zeros_like ( self.sig1_up ) # [kN/m]
        k_fl_NM_up = ones_like ( self.sig1_up ) # [-]
        k_fl_NM_lo = ones_like ( self.sig1_up ) # [-]
        alpha_up = zeros_like ( self.sig1_up ) # [deg]
        alpha_lo = zeros_like ( self.sig1_up ) # [deg]
        sig_comp_Ed_up = zeros_like ( self.sig1_up ) # [MPa]
        sig_comp_Ed_lo = zeros_like ( self.sig1_up ) # [MPa]

        #-------------------------------------------------
        # VAR 1:use simplified reinforced concrete approach
        #-------------------------------------------------
        if self.eval_mode == 'massivbau':

            m_Eds = zeros_like ( self.sig1_up ) # [kNm/m]
            e = zeros_like ( self.sig1_up ) # [m]

            zs = self.zs
            z = self.z

            #---------------------------------------------------------
            # tension upper side (sig1_up > 0):
            #---------------------------------------------------------

            # pure tension case:
            #
            bool_arr = cond_up_t
            alpha_up[ bool_arr ] = self.alpha_sig_up[ bool_arr ]
            m = abs( self.m_sig_up[ bool_arr ] )
            n = self.n_sig_up[ bool_arr ]
            # excentricity
            e[ bool_arr ] = abs( m / n )
            # in case of pure tension in the cross section:
            f_t_sig_up[ bool_arr ] = n * ( zs[ bool_arr ] + e[ bool_arr ] ) / ( zs[ bool_arr ] + zs[ bool_arr ] )

            ### bending with tension at the upper side 
            #
            bool_arr = cond_up_tb
            alpha_up[ bool_arr ] = self.alpha_sig_up[ bool_arr ]
            m = abs( self.m_sig_up[ bool_arr ] )
            n = self.n_sig_up[ bool_arr ]
            # moment at the height of the resulting reinforcement layer:
            m_Eds[ bool_arr ] = abs( m ) - zs[ bool_arr ] * n
            # tensile force in the reinforcement for bending and compression
            f_t_sig_up[ bool_arr ] = m_Eds[ bool_arr ] / z[ bool_arr ] + n

            #---------------------------------------------------------
            # compression upper side case (sig1_up < 0):
            #---------------------------------------------------------

            # @todo: remove as this is redundant:
            #
            bool_arr = cond_up_c
            alpha_up[ bool_arr ] = self.alpha_sig_up[ bool_arr ]
            f_t_sig_up[ bool_arr ] = 0.

            bool_arr = cond_up_cb
            alpha_up[ bool_arr ] = self.alpha_sig_up[ bool_arr ]
            f_t_sig_up[ bool_arr ] = 0.

            #---------------------------------------------------------
            # tension lower side (sig1_lo > 0):
            #---------------------------------------------------------

            # pure tension case:
            #
            bool_arr = cond_lo_t
            alpha_lo[ bool_arr ] = self.alpha_sig_lo[ bool_arr ]
            m = abs( self.m_sig_lo[ bool_arr ] )
            n = self.n_sig_lo[ bool_arr ]
            # excentricity
            e[ bool_arr ] = abs( m / n )
            # in case of pure tension in the cross section:
            f_t_sig_lo[ bool_arr ] = n * ( zs[ bool_arr ] + e[ bool_arr ] ) / ( zs[ bool_arr ] + zs[ bool_arr ] )

            ### bending with tension at the lower side 
            #
            bool_arr = cond_lo_tb
            alpha_lo[ bool_arr ] = self.alpha_sig_lo[ bool_arr ]
            m = abs( self.m_sig_lo[ bool_arr ] )
            n = self.n_sig_lo[ bool_arr ]
            # moment at the height of the resulting reinforcement layer:
            m_Eds[ bool_arr ] = abs( m ) - zs[ bool_arr ] * n
            # tensile force in the reinforcement for bending and compression
            f_t_sig_lo[ bool_arr ] = m_Eds[ bool_arr ] / z[ bool_arr ] + n

            #---------------------------------------------------------
            # compression lower side case (sig1_lo < 0):
            #---------------------------------------------------------

            # @todo: remove as this is redundant:
            #
            bool_arr = cond_lo_c
            alpha_lo[ bool_arr ] = self.alpha_sig_lo[ bool_arr ]
            f_t_sig_lo[ bool_arr ] = 0.

            bool_arr = cond_lo_cb
            alpha_lo[ bool_arr ] = self.alpha_sig_lo[ bool_arr ]
            f_t_sig_lo[ bool_arr ] = 0.


        #-------------------------------------------------
        # VAR 2:use principal stresses to calculate the resulting tensile force
        #-------------------------------------------------
        #
        if self.eval_mode == 'princ_sig_level_1':

            # NOTE: this evaluation mode has been used for the evaluation of 'n_tex' 
            # as described in the 'ZiE'. Note that in the evaluation of the principle 
            # stresses at the upper side a wrong variable has been used (had a minor effect) 
            # 

            print "NOTE: the principle tensile stresses are used to evaluate 'n_tex'"
            # resulting tensile force of the composite cross section[kN]
            # the entire (!) cross section is used!
            # as the maximum value of the tensile stresses at the top or the bottom 
            # i.e. sig1_max = min( 0, max( self.sig1_up, self.sig1_lo ) )

            ###-----------------------------------------------------------------------------------------------------------
            ### dimensioning ###
            ###-----------------------------------------------------------------------------------------------------------

            #---------------------------------------------------------
            # tension upper side (sig1_up > 0):
            #---------------------------------------------------------

            # pure tension case:
            # (--> 'k_fl_NM' set to 1.0, i.e. no increase of resistance due to bending)
            #
            bool_arr = cond_up_t
            alpha_up[ bool_arr ] = self.alpha_sig_up[ bool_arr ]
            sig_comp_Ed_up[ bool_arr ] = self.sig1_up[ bool_arr ]
            k_fl_NM_up[ bool_arr ] = 1.0
            f_t_sig_up[ bool_arr ] = self.sig1_up[ bool_arr ] * self.thickness[ bool_arr ] * 1000.

            # bending case with tension:
            # (sig_N > 0, i.e. 'sig1_up' results from bending and tension --> interpolate 'k_fl_NM')
            #
            bool_arr = cond_up_tb_t
            alpha_up[ bool_arr ] = self.alpha_sig_up[ bool_arr ]
            sig_comp_Ed_up[ bool_arr ] = self.sig1_up[ bool_arr ]
            sig_b = ( abs( self.sig1_up ) + abs( self.sig1_lo_sig_up ) ) / 2
            k_fl_NM_up[ bool_arr ] = 1.0 + ( self.k_fl - 1.0 ) * \
                                 ( 1.0 - ( sig_b[ bool_arr ] - abs( self.sig1_lo_sig_up[ bool_arr] ) ) / self.sig1_up[ bool_arr ] )
            f_t_sig_up[ bool_arr ] = self.sig1_up[ bool_arr ] * self.thickness[ bool_arr ] * 1000.

            # bending case with compression:
            # (sig_N < 0, i.e. 'sig1_up' results only from bending --> full increase for 'k_fl_NM')
            #
            bool_arr = cond_up_tb_c
            alpha_up[ bool_arr ] = self.alpha_sig_up[ bool_arr ]
            sig_comp_Ed_up[ bool_arr ] = self.sig1_up[ bool_arr ]
            sig_b = ( abs( self.sig1_up ) + abs( self.sig1_lo_sig_up ) ) / 2
            k_fl_NM_up[ bool_arr ] = self.k_fl
            f_t_sig_up[ bool_arr ] = self.sig1_up[ bool_arr ] * self.thickness[ bool_arr ] * 1000.

            #---------------------------------------------------------
            # compression upper side case (sig1_up < 0):
            #---------------------------------------------------------

            # @todo: remove as this is redundant:
            #
            bool_arr = cond_up_c
            alpha_up[ bool_arr ] = self.alpha_sig_up[ bool_arr ]
            sig_comp_Ed_up[ bool_arr ] = 0.
            k_fl_NM_up[ bool_arr ] = 1.
            f_t_sig_up[ bool_arr ] = 0.

            bool_arr = cond_up_cb
            alpha_up[ bool_arr ] = self.alpha_sig_up[ bool_arr ]
            sig_comp_Ed_up[ bool_arr ] = 0.
            k_fl_NM_up[ bool_arr ] = 1.
            f_t_sig_up[ bool_arr ] = 0.


            #---------------------------------------------------------
            # tension lower side (sig1_lo > 0):
            #---------------------------------------------------------

            # pure tension case:
            # (--> 'k_fl_NM' set to 1.0, i.e. no increase of resistance due to bending)
            #
            bool_arr = cond_lo_t
            alpha_lo[ bool_arr ] = self.alpha_sig_lo[ bool_arr ]
            sig_comp_Ed_lo[ bool_arr ] = self.sig1_lo[ bool_arr ]
            k_fl_NM_lo[ bool_arr ] = 1.0
            f_t_sig_lo[ bool_arr ] = self.sig1_lo[ bool_arr ] * self.thickness[ bool_arr ] * 1000.

            # bending case with tension:
            # (sig_N > 0, i.e. 'sig1_lo' results from bending and tension --> interpolate 'k_fl_NM')
            #
            bool_arr = cond_lo_tb_t
            alpha_lo[ bool_arr ] = self.alpha_sig_lo[ bool_arr ]
            sig_comp_Ed_lo[ bool_arr ] = self.sig1_lo[ bool_arr ]
            sig_b = ( abs( self.sig1_lo ) + abs( self.sig1_up_sig_lo ) ) / 2
            k_fl_NM_lo[ bool_arr ] = 1.0 + ( self.k_fl - 1.0 ) * \
                                 ( 1.0 - ( sig_b[ bool_arr ] - abs( self.sig1_up_sig_lo[ bool_arr] ) ) / self.sig1_lo[ bool_arr ] )
            f_t_sig_lo[ bool_arr ] = self.sig1_lo[ bool_arr ] * self.thickness[ bool_arr ] * 1000.

            # bending case with compression:
            # (sig_N < 0, i.e. 'sig1_lo' results only from bending --> full increase for 'k_fl_NM')
            #
            bool_arr = cond_lo_tb_c
            alpha_lo[ bool_arr ] = self.alpha_sig_lo[ bool_arr ]
            sig_comp_Ed_lo[ bool_arr ] = self.sig1_lo[ bool_arr ]
            sig_b = ( abs( self.sig1_lo ) + abs( self.sig1_up_sig_lo ) ) / 2
            k_fl_NM_lo[ bool_arr ] = self.k_fl
            f_t_sig_lo[ bool_arr ] = self.sig1_lo[ bool_arr ] * self.thickness[ bool_arr ] * 1000.

            #---------------------------------------------------------
            # compression lower side case (sig1_lo < 0):
            #---------------------------------------------------------

            # @todo: remove as this is redundant:
            #
            bool_arr = cond_lo_c
            alpha_lo[ bool_arr ] = self.alpha_sig_lo[ bool_arr ]
            sig_comp_Ed_lo[ bool_arr ] = 0.
            k_fl_NM_lo[ bool_arr ] = 1.
            f_t_sig_lo[ bool_arr ] = 0.

            bool_arr = cond_lo_cb
            alpha_lo[ bool_arr ] = self.alpha_sig_lo[ bool_arr ]
            sig_comp_Ed_lo[ bool_arr ] = 0.
            k_fl_NM_lo[ bool_arr ] = 1.
            f_t_sig_lo[ bool_arr ] = 0.


        #-------------------------------------------------
        # VAR 3:use principal stresses to calculate the resulting tensile force
        # using an internal iteration for the distribution of stresses
        #-------------------------------------------------
        #
        if self.eval_mode == 'princ_sig_level_2':

            print "NOTE: the principle tensile stresses are used to evaluate 'n_tex'"
            # resulting tensile force of the composite cross section[kN]
            # the entire (!) cross section is used!
            # as the maximum value of the tensile stresses at the top or the bottom 
            # i.e. sig1_max = min( 0, max( self.sig1_up, self.sig1_lo ) )

            #---------------------------------------------------------
            # prepare methods for bending cases to solve iteration 
            # to derive the tensile stresses of the outer layer
            #---------------------------------------------------------
            #
            sig_fl_calib = SigFlCalib( 
                               # concrete strength after 9 days
                               f_ck = 60.,
                               E_yarn = 84000,
                               width = 1.0,
                               n_roving = 120,
                               calib_config = 'eps_t_eps_c'
                              )

            def fit_response( m, n, thickness_unreinf, eps_t_estimate, eps_c_estimate ):
                u0_start = [eps_t_estimate, eps_c_estimate]
                return sig_fl_calib.fit_response( m, n, thickness_unreinf, u0_start )

            fit_response_vectorized = frompyfunc( fit_response, 5, 1 )

            get_sig_max_vectorized = frompyfunc( sig_fl_calib.get_sig_max, 2, 1 )

            # young's modulus of the composite cross section
            # used to derive a starting estimate of the composite stresses
            # in order to find converging solution in 'fit_response' 
            # @todo: derive that directly from 'sig_fl_calib'
            #
            E_comp = 30000.
            thickness_at_roof_edge = 0.060 #[m]

            #---------------------------------------------------------
            # tension upper side (sig1_up > 0):
            #---------------------------------------------------------

            # pure tension case:
            #
            bool_arr = cond_up_t
            alpha_up[ bool_arr ] = self.alpha_sig_up[ bool_arr ]
            sig_comp_Ed_up[ bool_arr ] = self.sig1_up[ bool_arr ]
            f_t_sig_up[ bool_arr ] = self.sig1_up[ bool_arr ] * self.thickness[ bool_arr ] * 1000.

            # bending case:
            #
            bool_arr = cond_up_tb
            alpha_up[ bool_arr ] = self.alpha_sig_up[ bool_arr ]
            #    
            m = abs( self.m_sig_up[ bool_arr ] )
            n = self.n_sig_up[ bool_arr ]
            #
            thickness_unreinf = self.thickness[ bool_arr ] - thickness_at_roof_edge
            eps_t_estimate = abs( self.sig1_up[ bool_arr ] ) / E_comp
            eps_c_estimate = abs( self.sig1_lo_sig_up[ bool_arr ] ) / E_comp
#            eps_t_estimate = 0.003
#            eps_c_estimate = 0.010
            u_sol_arr = fit_response_vectorized( m, n, thickness_unreinf, eps_t_estimate, eps_c_estimate )
            max_sig_arr = get_sig_max_vectorized( u_sol_arr, thickness_unreinf )
            #
            sig_comp_Ed_up[ bool_arr ] = max_sig_arr
            f_t_sig_up[ bool_arr ] = self.sig1_up[ bool_arr ] * self.thickness[ bool_arr ] * 1000.

            #---------------------------------------------------------
            # compression upper side case (sig1_up < 0):
            #---------------------------------------------------------

            # @todo: remove as this is redundant:
            #
            bool_arr = cond_up_c
            alpha_up[ bool_arr ] = self.alpha_sig_up[ bool_arr ]
            sig_comp_Ed_up[ bool_arr ] = 0.
            f_t_sig_up[ bool_arr ] = 0.

            bool_arr = cond_up_cb
            alpha_up[ bool_arr ] = self.alpha_sig_up[ bool_arr ]
            sig_comp_Ed_up[ bool_arr ] = 0.
            f_t_sig_up[ bool_arr ] = 0.


            #---------------------------------------------------------
            # tension lower side (sig1_lo > 0):
            #---------------------------------------------------------

            # pure tension case:
            #
            bool_arr = cond_lo_t
            alpha_lo[ bool_arr ] = self.alpha_sig_lo[ bool_arr ]
            sig_comp_Ed_lo[ bool_arr ] = self.sig1_lo[ bool_arr ]
            f_t_sig_lo[ bool_arr ] = self.sig1_lo[ bool_arr ] * self.thickness[ bool_arr ] * 1000.

            # bending case:
            #
            bool_arr = cond_lo_tb
            alpha_lo[ bool_arr ] = self.alpha_sig_lo[ bool_arr ]
            #    
            m = abs( self.m_sig_lo[ bool_arr ] )
            n = self.n_sig_lo[ bool_arr ]
            #
            thickness_unreinf = self.thickness[ bool_arr ] - thickness_at_roof_edge
            eps_t_estimate = abs( self.sig1_lo[ bool_arr ] ) / E_comp
            eps_c_estimate = abs( self.sig1_up_sig_lo[ bool_arr ] ) / E_comp
#            eps_t_estimate = 0.003
#            eps_c_estimate = 0.010
            u_sol_arr = fit_response_vectorized( m, n, thickness_unreinf, eps_t_estimate, eps_c_estimate )
            max_sig_arr = get_sig_max_vectorized( u_sol_arr, thickness_unreinf )
            #
            sig_comp_Ed_lo[ bool_arr ] = max_sig_arr.flatten()
            f_t_sig_lo[ bool_arr ] = self.sig1_lo[ bool_arr ] * self.thickness[ bool_arr ] * 1000.

            #---------------------------------------------------------
            # compression lower side case (sig1_lo < 0):
            #---------------------------------------------------------

            # @todo: remove as this is redundant:
            #
            bool_arr = cond_lo_c
            alpha_lo[ bool_arr ] = self.alpha_sig_lo[ bool_arr ]
            sig_comp_Ed_lo[ bool_arr ] = 0.
            f_t_sig_lo[ bool_arr ] = 0.

            bool_arr = cond_lo_cb
            alpha_lo[ bool_arr ] = self.alpha_sig_lo[ bool_arr ]
            sig_comp_Ed_lo[ bool_arr ] = 0.
            f_t_sig_lo[ bool_arr ] = 0.


            #---------------------------------------------------------


        #------------------------------------------------------------
        # get angel of deflection of the textile reinforcement
        #------------------------------------------------------------
        # angel of deflection of the textile reinforcement 
        # distinguished between longitudinal (l) and transversal (q) direction
        print "NOTE: deflection angle is used to evaluate 'n_tex'"

        beta_l_up_deg = abs( alpha_up ) # [degree]   
        beta_q_up_deg = 90. - abs( alpha_up )      # [degree]   
        beta_l_up = beta_l_up_deg * pi / 180. # [rad] 
        beta_q_up = beta_q_up_deg * pi / 180. # [rad]

        beta_l_lo_deg = 90 - abs( alpha_lo ) # [degree]   
        beta_q_lo_deg = abs( alpha_lo )      # [degree]   
        beta_l_lo = beta_l_lo_deg * pi / 180. # [rad] 
        beta_q_lo = beta_q_lo_deg * pi / 180. # [rad]

        # resulting strength of the bi-directional textile considering the 
        # deflection of the reinforcement in the loading direction
        # per reinforcement layer
        #
        f_Rtex_0 = self.f_Rtex_0  # [kN/m/layer]
        f_Rtex_90 = self.f_Rtex_90
        f_Rtex_lo = f_Rtex_0 * cos( beta_l_lo ) * ( 1 - beta_l_lo_deg / 90. ) + \
                    f_Rtex_90 * cos( beta_q_lo ) * ( 1 - beta_q_lo_deg / 90. )
        f_Rtex_up = f_Rtex_0 * cos( beta_l_up ) * ( 1 - beta_l_up_deg / 90. ) + \
                    f_Rtex_90 * cos( beta_q_up ) * ( 1 - beta_q_up_deg / 90. )

        # @todo: check for general case 
        # NOTE: needs information about the orientation of the reinforcement
        # works here only because of the simplification that the same resistance of the textile in 0- and 90-direction is assumed
        # and the reinforcement is arranged in the shell approximately orthogonal to
        # the global coordinate system
        #
        sig_comp_Rd_lo = self.sig_comp_0_Rd * cos( beta_l_lo ) * ( 1 - beta_l_lo_deg / 90. ) + \
                         self.sig_comp_90_Rd * cos( beta_q_lo ) * ( 1 - beta_q_lo_deg / 90. )
        sig_comp_Rd_up = self.sig_comp_0_Rd * cos( beta_l_up ) * ( 1 - beta_l_up_deg / 90. ) + \
                         self.sig_comp_90_Rd * cos( beta_q_up ) * ( 1 - beta_q_up_deg / 90. )


        #------------------------------------------------------------
        # construct a dictionary containing the return values
        #------------------------------------------------------------

        if self.eval_mode == 'massivbau':

            # necessary number of reinforcement layers
            # for the entire cross section (use symmetric arrangement of the upper and
            # lower reinforcement layers
            #
            n_tex_up = f_t_sig_up / f_Rtex_up
            n_tex_lo = f_t_sig_lo / f_Rtex_lo
            n_tex = 2 * ndmax( hstack( [ n_tex_up, n_tex_lo] ), axis = 1 )[:, None]

            return { 'e':e, 'm_Eds':m_Eds,
                     'cond_up_tb' : cond_up_tb * 1.0,
                     'cond_lo_tb' : cond_lo_tb * 1.0,
                     'f_t_sig_up' : f_t_sig_up,
                     'f_t_sig_lo' : f_t_sig_lo,
                     'beta_l_up':beta_l_up_deg, 'beta_q_up':beta_q_up_deg,
                     'beta_l_lo':beta_l_lo_deg, 'beta_q_lo':beta_q_lo_deg,
                     'f_Rtex_up':f_Rtex_up,
                     'f_Rtex_lo':f_Rtex_lo,
                     'n_tex_up':n_tex_up,
                     'n_tex_lo':n_tex_lo,
                     'n_tex':n_tex}



        if self.eval_mode == 'princ_sig_level_1':

            # ratio of the imposed stresses and the composite resistance
            # NOTE: resistance is increased by factor 'k_fl_NM' if a bending case is evaluated
            #
            eta_comp_up = sig_comp_Ed_up / ( sig_comp_Rd_up * k_fl_NM_up )
            eta_comp_lo = sig_comp_Ed_lo / ( sig_comp_Rd_lo * k_fl_NM_lo )
            eta_comp = ndmax( hstack( [ eta_comp_up, eta_comp_lo] ), axis = 1 )[:, None]

            # necessary number of reinforcement layers
            # for the cases that stresses at the upper or lower face are taken into account
            # for the evaluation of the necessary number of reinforcement layers. 'n_tex' is the
            # maximum of the upper and lower face evaluation. 
            # NOTE: resistance is increased by factor 'k_fl_NM' if a bending case is evaluated
            #
            n_tex_up = f_t_sig_up / ( f_Rtex_up * k_fl_NM_up )
            n_tex_lo = f_t_sig_lo / ( f_Rtex_lo * k_fl_NM_lo )

            # use a symmetric reinforcement layup at the top and at the bottom:
            #
            n_tex = ndmax( hstack( [ n_tex_up, n_tex_lo] ), axis = 1 )[:, None]

            return {
                     'cond_up_tb' : cond_up_tb * 1.0,
                     'cond_lo_tb' : cond_lo_tb * 1.0,
                     'f_t_sig_up' : f_t_sig_up,
                     'f_t_sig_lo' : f_t_sig_lo,
                     'beta_l_up':beta_l_up_deg, 'beta_q_up':beta_q_up_deg,
                     'beta_l_lo':beta_l_lo_deg, 'beta_q_lo':beta_q_lo_deg,
                     'sig_comp_Rd_up':sig_comp_Rd_up,
                     'sig_comp_Rd_lo':sig_comp_Rd_lo,
                     'eta_comp_up':eta_comp_up,
                     'eta_comp_lo':eta_comp_lo,
                     'eta_comp':eta_comp,
                     'f_Rtex_up':f_Rtex_up,
                     'k_fl_NM_up':k_fl_NM_up,
                     'f_Rtex_lo':f_Rtex_lo,
                     'k_fl_NM_lo':k_fl_NM_lo,
                     'n_tex_up':n_tex_up,
                     'n_tex_lo':n_tex_lo,
                     'n_tex':n_tex}

        if self.eval_mode == 'princ_sig_level_2':

            # ratio of the imposed stresses and the composite resistance
            #
            eta_comp_up = sig_comp_Ed_up / sig_comp_Rd_up
            eta_comp_lo = sig_comp_Ed_lo / sig_comp_Rd_lo
            eta_comp = ndmax( hstack( [ eta_comp_up, eta_comp_lo] ), axis = 1 )[:, None]

            n_tex_up = f_t_sig_up / f_Rtex_up
            n_tex_lo = f_t_sig_lo / f_Rtex_lo

            # use a symmetric reinforcement layup at the top and at the bottom:
            #
            n_tex = ndmax( hstack( [ n_tex_up, n_tex_lo] ), axis = 1 )[:, None]

            return {
                     'cond_up_tb' : cond_up_tb * 1.0,
                     'cond_lo_tb' : cond_lo_tb * 1.0,
                     'f_t_sig_up' : f_t_sig_up,
                     'f_t_sig_lo' : f_t_sig_lo,
                     'beta_l_up':beta_l_up_deg, 'beta_q_up':beta_q_up_deg,
                     'beta_l_lo':beta_l_lo_deg, 'beta_q_lo':beta_q_lo_deg,
                     'sig_comp_Rd_up':sig_comp_Rd_up,
                     'sig_comp_Rd_lo':sig_comp_Rd_lo,
                     'eta_comp_up':eta_comp_up,
                     'eta_comp_lo':eta_comp_lo,
                     'eta_comp':eta_comp,
                     'f_Rtex_up':f_Rtex_up,
                     'f_Rtex_lo':f_Rtex_lo,
                     'n_tex_up':n_tex_up,
                     'n_tex_lo':n_tex_lo,
                     'n_tex':n_tex}


    #-----------------------------------------------
    # LS_COLUMNS: specify the properties that are displayed in the view
    #-----------------------------------------------

    if eval_mode == 'massivbau':

        assess_name = 'max_n_tex'

        ls_columns = List( ['d', 'zs', 'z',
                            'e', 'm_Eds',
                            'f_t_sig_up', 'f_t_sig_lo',
                            'beta_l_up', 'beta_q_up',
                            'beta_l_lo', 'beta_q_lo',
                            'f_Rtex_up', 'f_Rtex_lo',
                            'n_tex_up', 'n_tex_lo',
                            'n_tex' ] )

        e = Property( Array )
        def _get_e( self ):
            return self.ls_values['e']

        m_Eds = Property( Array )
        def _get_m_Eds( self ):
            return self.ls_values['m_Eds']

    elif eval_mode == 'princ_sig_level_1':

#        assess_name = 'max_eta_comp'
        assess_name = 'max_n_tex'

        ls_columns = List( ['f_t_sig_up', 'f_t_sig_lo',
                            'beta_l_up', 'beta_q_up',
                            'beta_l_lo', 'beta_q_lo',
                            'f_Rtex_up', 'f_Rtex_lo',
                            'k_fl_NM_up', 'k_fl_NM_lo',
                            'eta_comp_up', 'eta_comp_lo', 'eta_comp',
                            'cond_up_tb', 'cond_lo_tb',
                            'n_tex_up', 'n_tex_lo', 'n_tex'] )

    elif eval_mode == 'princ_sig_level_2':

        assess_name = 'max_eta_comp'

        ls_columns = List( ['f_t_sig_up', 'f_t_sig_lo',
                            'beta_l_up', 'beta_q_up',
                            'beta_l_lo', 'beta_q_lo',
                            'f_Rtex_up', 'f_Rtex_lo',
                            'eta_comp_up', 'eta_comp_lo', 'eta_comp',
                            'cond_up_tb', 'cond_lo_tb',
                            'n_tex_up', 'n_tex_lo', 'n_tex'] )


    k_fl_NM_lo = Property( Array )
    def _get_k_fl_NM_lo( self ):
        return self.ls_values['k_fl_NM_lo']

    k_fl_NM_up = Property( Array )
    def _get_k_fl_NM_up( self ):
        return self.ls_values['k_fl_NM_up']

    cond_up_tb = Property( Array )
    def _get_cond_up_tb( self ):
        return self.ls_values['cond_up_tb']

    cond_lo_tb = Property( Array )
    def _get_cond_lo_tb( self ):
        return self.ls_values['cond_lo_tb']

    f_t_sig_up = Property( Array )
    def _get_f_t_sig_up( self ):
        return self.ls_values['f_t_sig_up']

    f_t_sig_lo = Property( Array )
    def _get_f_t_sig_lo( self ):
        return self.ls_values['f_t_sig_lo']

    beta_l_up = Property( Array )
    def _get_beta_l_up( self ):
        return self.ls_values['beta_l_up']

    beta_q_up = Property( Array )
    def _get_beta_q_up( self ):
        return self.ls_values['beta_q_up']

    beta_l_lo = Property( Array )
    def _get_beta_l_lo( self ):
        return self.ls_values['beta_l_lo']

    beta_q_lo = Property( Array )
    def _get_beta_q_lo( self ):
        return self.ls_values['beta_q_lo']

    f_Rtex_up = Property( Array )
    def _get_f_Rtex_up( self ):
        return self.ls_values['f_Rtex_up']

    f_Rtex_lo = Property( Array )
    def _get_f_Rtex_lo( self ):
        return self.ls_values['f_Rtex_lo']

    k_fl_NM = Property( Array )
    def _get_k_fl_NM( self ):
        return self.ls_values['k_fl_NM']

    n_tex_up = Property( Array )
    def _get_n_tex_up( self ):
        return self.ls_values['n_tex_up']

    n_tex_lo = Property( Array )
    def _get_n_tex_lo( self ):
        return self.ls_values['n_tex_lo']

    n_tex = Property( Array )
    def _get_n_tex( self ):
        return self.ls_values['n_tex']

    eta_comp_up = Property( Array )
    def _get_eta_comp_up( self ):
        return self.ls_values['eta_comp_up']

    eta_comp_lo = Property( Array )
    def _get_eta_comp_lo( self ):
        return self.ls_values['eta_comp_lo']

    eta_comp = Property( Array )
    def _get_eta_comp( self ):
        return self.ls_values['eta_comp']


    #-------------------------------------------------
    # choose the assess parameter used for sorting
    # defined by the property name 'assess_name'
    #-------------------------------------------------

#    assess_name = 'max_eta_comp'
#    assess_name = 'max_n_tex'

    max_eta_comp = Property( depends_on = '+input' )
    @cached_property
    def _get_max_eta_comp( self ):
        return ndmax( self.eta_comp )

    max_n_tex = Property( depends_on = '+input' )
    @cached_property
    def _get_max_n_tex( self ):
        return ndmax( self.n_tex )

#    assess_name = 'max_sig1_up'
#    max_sig1_up = Property( depends_on = '+input' )
#    @cached_property
#    def _get_max_sig1_up( self ):
#        return ndmax( self.sig1_up )

#    assess_name = 'max_sig1_lo'
#    max_sig1_lo = Property( depends_on = '+input' )
#    @cached_property
#    def _get_max_sig1_lo( self ):
#        return ndmax( self.sig1_lo )


    #@todo: make an automatised function calle for asses_value_max
#    @on_trait_change( '+input' )
#    def set_assess_name_max( self, assess_name ):
#        print 'set asses'
#        asses_value = ndmax( getattr( self, assess_name ) )
#        assess_name_max = 'max_' + assess_name
#        setattr( self, assess_name_max, asses_value )




    #-------------------------------
    # ls view
    #-------------------------------

    # @todo: the dynamic selection of the columns to be displayed 
    # does not work in connection with the LSArrayAdapter 
    traits_view = View( 
                       VGroup( 
                        HGroup( 
                            VGroup( 
                                Item( name = 'sig_comp_0_Rd', label = 'composit tensile strength [MPa]:  sig_comp_0_Rd ', style = 'readonly', format_str = "%.1f" ),
                                Item( name = 'f_Rtex_0', label = 'reinforcement strength per layer [kN/m]:  f_Rtex_0 ', style = 'readonly', format_str = "%.1f" ),
                                label = 'material properties (longitudinal)'
                                  ),
                            VGroup( 
                                Item( name = 'sig_comp_90_Rd', label = 'composit tensile strength [MPa]:  sig_comp_90_Rd ', style = 'readonly', format_str = "%.1f" ),
                                Item( name = 'f_Rtex_90', label = 'reinforcement strength per layer [kN/m]:  f_Rtex_90 ', style = 'readonly', format_str = "%.1f" ),
                                label = 'material Properties (transversal)'
                                  ),
                             ),

                        VGroup( 
                            Include( 'ls_group' ),
                            Item( 'ls_array', show_label = False,
                                  editor = TabularEditor( adapter = LSArrayAdapter() ) )
                              ),
                            ),
                      resizable = True,
                      scrollable = True,
                      height = 1000,
                      width = 1100
                      )

LSLIST = [ SLS, ULS ]

class LSTable( HasTraits ):
    '''Assessment tool
    '''

    is_id = Int( 0 )
    # geo data: coordinates and element thickness
    # 
    geo_data = Dict

    elem_no = Property( Array )
    def _get_elem_no( self ):
        return self.geo_data['elem_no']

    X = Property( Array )
    def _get_X( self ):
        return self.geo_data['X']

    Y = Property( Array )
    def _get_Y( self ):
        return self.geo_data['Y']

    Z = Property( Array )
    def _get_Z( self ):
        return self.geo_data['Z']

    thickness = Property( Array )
    def _get_thickness( self ):
        '''element thickness (units changed form [mm] to [m])'''
        return self.geo_data['thickness'] / 1000.

    # state data: stress resultants 
    # 
    state_data = Dict

    mx = Property( Array )
    def _get_mx( self ):
        return self.state_data['mx']

    my = Property( Array )
    def _get_my( self ):
        return self.state_data['my']

    mxy = Property( Array )
    def _get_mxy( self ):
        return self.state_data['mxy']

    nx = Property( Array )
    def _get_nx( self ):
        return self.state_data['nx']

    ny = Property( Array )
    def _get_ny( self ):
        return self.state_data['ny']

    nxy = Property( Array )
    def _get_nxy( self ):
        return self.state_data['nxy']

    # ------------------------------------------------------------
    # Index sig: calculate principle direction of the stresses at 
    # the lower and upper side and get the corresponding values of
    # the stresses at the opposite side. Also get the corresponding
    # values of the normal force and the moment in this direction
    # ------------------------------------------------------------

    princ_values_sig = Property( Dict, depends_on = 'data_file_stress_resultants' )
    @cached_property
    def _get_princ_values_sig( self ):
        '''principle value of the stresses for the lower ('lo') and upper ('up') face:
        '''
        # stress_resultants in global coordinates
        # --> moments in kNm
        # --> normal forces in kN
        # --> area in m**2
        # --> resisting moment in m**3
        # --> stresses in MPa
        #
        mx = self.mx
        my = self.my
        mxy = self.mxy
        nx = self.nx
        ny = self.ny
        nxy = self.nxy

        # geometrical properties:
        #
        A = self.thickness * 1.0
        W = self.thickness ** 2 * 1.0 / 6.


        # compare the formulae with the RFEM-manual p.290

        # stresses [MPa] upper face in global drection:
        #
        sigx_up = ( nx / A - mx / W ) / 1000.
        sigy_up = ( ny / A - my / W ) / 1000.
        sigxy_up = ( nxy / A - mxy / W ) / 1000.

        # stresses [MPa] lower face in global direction:
        #
        sigx_lo = ( nx / A + mx / W ) / 1000.
        sigy_lo = ( ny / A + my / W ) / 1000.
        sigxy_lo = ( nxy / A + mxy / W ) / 1000.


        # principal stresses upper face:
        #
        sig1_up = 0.5 * ( sigx_up + sigy_up ) + 0.5 * sqrt( ( sigx_up - sigy_up ) ** 2 + 4 * sigxy_up ** 2 )
        sig2_up = 0.5 * ( sigx_up + sigy_up ) - 0.5 * sqrt( ( sigx_up - sigy_up ) ** 2 + 4 * sigxy_up ** 2 )

        alpha_sig_up = pi / 2. * ones_like( sig1_up )

        # from mechanic formula book (cf. also InfoCAD manual)
        bool = sig2_up != sigx_up

        ### NOTE: for ZiE the following line was used containing one wrong parameter ('sigy_up' instead of 'sig2_up')
        ### this leads to conservative values (more reinforcement then necessary  --> n_tex = 9.78 insted of 9.42)
        #
#        alpha_sig_up[ bool ] = arctan( sigxy_up[ bool ] / ( sigy_up[ bool ] - sigx_up[ bool ] ) )

        alpha_sig_up[ bool ] = arctan( sigxy_up[ bool ] / ( sig2_up[ bool ] - sigx_up[ bool ] ) )

        # RFEM-manual (NOTE that manual contains typing error!)
        # the formula as given below yields the same results then the used mechanic formula
#        bool = sigx_up != sigy_up
#        alpha_sig_up[ bool ] = 0.5 * arctan( 2 * sigxy_up[ bool ] / ( sigx_up[ bool ] - sigy_up[ bool ] ) )

        alpha_sig_up_deg = alpha_sig_up * 180. / pi

        # transform formula taken from mechanic formula book
        # transform the stresses at the lower face to the principle tensile direction (1-direction) of the upper stresses
        sig1_lo_sig_up = 0.5 * ( sigy_lo + sigx_lo ) - 0.5 * ( sigy_lo - sigx_lo ) * cos( 2 * alpha_sig_up ) - sigxy_lo * sin( 2 * alpha_sig_up )

        # transform moments and normal forces in the direction of the principal stresses (1-direction)
        #
        m_sig_up = 0.5 * ( my + mx ) - 0.5 * ( my - mx ) * cos( 2 * alpha_sig_up ) - mxy * sin( 2 * alpha_sig_up )
        n_sig_up = 0.5 * ( ny + nx ) - 0.5 * ( ny - nx ) * cos( 2 * alpha_sig_up ) - nxy * sin( 2 * alpha_sig_up )


        #--------------
        # lower face:
        #--------------

        # principal stresses lower face:
        #
        sig1_lo = 0.5 * ( sigx_lo + sigy_lo ) + 0.5 * sqrt( ( sigx_lo - sigy_lo ) ** 2 + 4 * sigxy_lo ** 2 )
        sig2_lo = 0.5 * ( sigx_lo + sigy_lo ) - 0.5 * sqrt( ( sigx_lo - sigy_lo ) ** 2 + 4 * sigxy_lo ** 2 )

        alpha_sig_lo = pi / 2. * ones_like( sig1_lo )

        # from mechanic formula book (cf. also InfoCAD manual)
        bool = sig2_lo != sigx_lo
        alpha_sig_lo[ bool ] = arctan( sigxy_lo[ bool ] / ( sig2_lo[ bool ] - sigx_lo[ bool ] ) )

        # RFEM-manual (NOTE that manual contains typing error!) 
        # the formula as given below yields the same results then the used mechanic formula
#        bool = sigx_lo != sigy_lo
#        alpha_sig_lo[ bool ] = 0.5 * arctan( 2 * sigxy_lo[ bool ] / ( sigx_lo[ bool ] - sigy_lo[ bool ] ) )

        alpha_sig_lo_deg = alpha_sig_lo * 180. / pi

        # transform the stresses at the lower face to the principle tensile direction (1-direction) of the upper stresses
        # Note: transformation forumla taken from mechanic formula book
        #
        sig1_up_sig_lo = 0.5 * ( sigy_up + sigx_up ) - 0.5 * ( sigy_up - sigx_up ) * cos( 2 * alpha_sig_lo ) - sigxy_up * sin( 2 * alpha_sig_lo )

        # transform moments and normal forces in the direction of the principal stresses (1-direction)
        #
        m_sig_lo = 0.5 * ( my + mx ) - 0.5 * ( my - mx ) * cos( 2 * alpha_sig_lo ) - mxy * sin( 2 * alpha_sig_lo )
        n_sig_lo = 0.5 * ( ny + nx ) - 0.5 * ( ny - nx ) * cos( 2 * alpha_sig_lo ) - nxy * sin( 2 * alpha_sig_lo )

        return {
                 'sigx_up' : sigx_up, 'sigy_up' : sigy_up, 'sigxy_up' : sigxy_up,
                 'sig1_up' : sig1_up, 'sig2_up' : sig2_up, 'alpha_sig_up' : alpha_sig_up_deg,

                 'sig1_lo_sig_up' : sig1_lo_sig_up,
                 'm_sig_up' : m_sig_up, 'n_sig_up' : n_sig_up,

                 'sigx_lo' : sigx_lo, 'sigy_lo' : sigy_lo, 'sigxy_lo' : sigxy_lo,
                 'sig1_lo' : sig1_lo, 'sig2_lo' : sig2_lo, 'alpha_sig_lo' : alpha_sig_lo_deg,

                 'sig1_up_sig_lo' : sig1_up_sig_lo,
                 'm_sig_lo' : m_sig_lo, 'n_sig_lo' : n_sig_lo,

               }

    # stresses upper face:
    #
    sigx_up = Property( Float )
    def _get_sigx_up( self ):
        return self.princ_values_sig['sigx_up']

    sigy_up = Property( Float )
    def _get_sigy_up( self ):
        return self.princ_values_sig['sigy_up']

    sigxy_up = Property( Float )
    def _get_sigxy_up( self ):
        return self.princ_values_sig['sigxy_up']

    sig1_up = Property( Float )
    def _get_sig1_up( self ):
        return self.princ_values_sig['sig1_up']

    sig2_up = Property( Float )
    def _get_sig2_up( self ):
        return self.princ_values_sig['sig2_up']

    alpha_sig_up = Property( Float )
    def _get_alpha_sig_up( self ):
        return self.princ_values_sig['alpha_sig_up']

    sig1_lo_sig_up = Property( Float )
    def _get_sig1_lo_sig_up( self ):
        return self.princ_values_sig['sig1_lo_sig_up']

    m_sig_up = Property( Float )
    def _get_m_sig_up( self ):
        return self.princ_values_sig['m_sig_up']

    n_sig_up = Property( Float )
    def _get_n_sig_up( self ):
        return self.princ_values_sig['n_sig_up']


    # stresses lower face:
    #
    sigx_lo = Property( Float )
    def _get_sigx_lo( self ):
        return self.princ_values_sig['sigx_lo']

    sigy_lo = Property( Float )
    def _get_sigy_lo( self ):
        return self.princ_values_sig['sigy_lo']

    sigxy_lo = Property( Float )
    def _get_sigxy_lo( self ):
        return self.princ_values_sig['sigxy_lo']

    sig1_lo = Property( Float )
    def _get_sig1_lo( self ):
        return self.princ_values_sig['sig1_lo']

    sig2_lo = Property( Float )
    def _get_sig2_lo( self ):
        return self.princ_values_sig['sig2_lo']

    alpha_sig_lo = Property( Float )
    def _get_alpha_sig_lo( self ):
        return self.princ_values_sig['alpha_sig_lo']

    sig1_up_sig_lo = Property( Float )
    def _get_sig1_up_sig_lo( self ):
        return self.princ_values_sig['sig1_up_sig_lo']

    m_sig_lo = Property( Float )
    def _get_m_sig_lo( self ):
        return self.princ_values_sig['m_sig_lo']

    n_sig_lo = Property( Float )
    def _get_n_sig_lo( self ):
        return self.princ_values_sig['n_sig_lo']


    #------------------------------------------
    # combinations of limit states, stress resultants and directions
    #------------------------------------------

    ls = Trait( 'ULS',
                {'ULS' : ULS,
                 'SLS' : SLS } )

    ls_class = Instance( LS )
    def _ls_class_default( self ):
        '''ls instances, e.g. ULS()
        '''
        ls_class = self.ls_
        return ls_class( ls_table = self )

    assess_value = Property
    def _get_assess_value( self ):
        ls = self.ls_class
        return getattr( ls, ls.assess_name )

    traits_view = View( Tabbed( 
                            Item( 'ls_class@' , label = "ls", show_label = False ),
                            scrollable = False,
                         ),
                      resizable = True,
                      scrollable = True,
                      height = 1000,
                      width = 1100
                      )


