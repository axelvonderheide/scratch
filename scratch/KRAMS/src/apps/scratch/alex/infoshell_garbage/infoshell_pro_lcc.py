'''
Created on Jun 23, 2010

@author: alexander
'''

from enthought.traits.api import \
    HasTraits, Directory, List, Int, Float, Any, Enum, \
    on_trait_change, File, Constant, Instance, Trait, \
    Array, Str, Property, cached_property, WeakRef, \
    Dict, Button, Color, Bool

from enthought.util.home_directory import \
    get_home_directory

from enthought.traits.ui.api import \
    View, Item, DirectoryEditor, TabularEditor, HSplit, Tabbed, VGroup, \
    TableEditor, Group, ListEditor, VSplit, HSplit, VGroup, HGroup, Spring, \
    Include

from enthought.mayavi import \
    mlab

from enthought.traits.ui.table_column import \
    ObjectColumn

from enthought.traits.ui.menu import \
    OKButton, CancelButton

from enthought.traits.ui.tabular_adapter \
    import TabularAdapter

from numpy import array, loadtxt, arange, sqrt, zeros, arctan, sin, cos, ones_like, \
                vstack, savetxt, hstack, argsort, fromstring, zeros_like, \
                copy, c_, newaxis, argmax, where, argsort, sqrt, frompyfunc

from math import pi
from string import split
import os

from scipy.io import read_array

from lc_manager import \
    LC, LCManager

DIRLIST = ['x', 'y']
SRLIST = ['M', 'N']

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
    info_shell = WeakRef

    # parameters of the limit state
    #
    dir = Enum( DIRLIST )
    stress_res = Enum( SRLIST )

    #-------------------------------
    # ls columns
    #-------------------------------
    # defined in the subclasses
    #
    ls_columns = List
    show_ls_columns = Bool( True )

    #-------------------------------
    # sr columns
    #-------------------------------

    # stress resultant columns - for ULS this is defined in the subclasses
    #
    sr_columns = List( ['m', 'n'] )
    show_sr_columns = Bool( True )

    # stress resultant columns - generated from the parameter combination
    # dir and stress_res - one of 'Mx', 'Ny', 'Mx', 'Ny'
    #
    m_varname = Property( Str )
    def _get_m_varname( self ):
        # e.g. mx_N 
        appendix = self.dir + '_' + self.stress_res
        return 'm' + appendix

    n_varname = Property( Str )
    def _get_n_varname( self ):
        # e.g. nx_N 
        appendix = self.dir + '_' + self.stress_res
        return 'n' + appendix

    n = Property( Float )
    def _get_n( self ):
        return getattr( self.info_shell, self.n_varname )

    m = Property( Float )
    def _get_m( self ):
        return getattr( self.info_shell, self.m_varname )

    #-------------------------------
    # geo columns form info shell
    #-------------------------------

    geo_columns = List( [ 'elem_no', 'X', 'Y', 'Z', 'D_elem' ] )
    show_geo_columns = Bool( True )

    elem_no = Property( Float )
    def _get_elem_no( self ):
        return self.info_shell.elem_no

    X = Property( Float )
    def _get_X( self ):
        return self.info_shell.X

    Y = Property( Float )
    def _get_Y( self ):
        return self.info_shell.Y

    Z = Property( Float )
    def _get_Z( self ):
        return self.info_shell.Z

    D_elem = Property( Float )
    def _get_D_elem( self ):
        return self.info_shell.D_elem

    #-------------------------------
    # state columns form info shell
    #-------------------------------

    state_columns = List( ['mx', 'my', 'mxy', 'nx', 'ny', 'nxy', 'combi_key' ] )
    show_state_columns = Bool( True )

    mx = Property( Float )
    def _get_mx( self ):
        return self.info_shell.mx

    my = Property( Float )
    def _get_my( self ):
        return self.info_shell.my

    mxy = Property( Float )
    def _get_mxy( self ):
        return self.info_shell.mxy

    nx = Property( Float )
    def _get_nx( self ):
        return self.info_shell.nx

    ny = Property( Float )
    def _get_ny( self ):
        return self.info_shell.ny

    nxy = Property( Float )
    def _get_nxy( self ):
        return self.info_shell.nxy

    combi_key = Property( Float )
    def _get_combi_key( self ):
        return self.info_shell.combi_key

    #-------------------------------
    # ls table
    #-------------------------------

    # all columns associated with the limit state including the corresponding
    # stress resultants
    #
    columns = Property( List, depends_on = 'show_geo_columns, show_state_columns,\
                                             show_sr_columns, show_ls_columns' )
    @cached_property
    def _get_columns( self ):
        columns = []

        if self.show_geo_columns:
            columns += self.geo_columns

        if self.show_state_columns:
            columns += self.state_columns

        if self.show_sr_columns:
            columns += self.sr_columns

        if self.show_ls_columns:
            columns += self.ls_columns

        return columns

    # select column used for sorting the data in selected sorting order 
    #
    sort_column = Enum( values = 'columns' )
    def _sort_column_default( self ):
        return self.columns[-1]

    sort_order = Enum( 'descending', 'ascending', 'unsorted' )

    # get the maximum value of the chosen column
    #
    max_in_column = Enum( values = 'columns' )
    def _max_in_column_default( self ):
        return self.columns[-1]

    max_value = Property( depends_on = 'max_in_column' )
    def _get_max_value( self ):
        col = getattr( self, self.max_in_column )[:, 0]
        return max( col )

    # stack columns together for table used by TabularEditor
    #
    ls_table = Property( Array, depends_on = 'sort_column, sort_order, show_geo_columns, \
                                              show_state_columns, show_sr_columns, show_ls_columns' )
    @cached_property
    def _get_ls_table( self ):

        arr_list = [ getattr( self, col ) for col in self.columns ]

        # get the array currently selected by the sort_column enumeration
        #
        sort_arr = getattr( self, self.sort_column )[:, 0]
        sort_idx = argsort( sort_arr )
        ls_table = hstack( arr_list )

        if self.sort_order == 'descending':
            return ls_table[ sort_idx[::-1] ]
        if self.sort_order == 'ascending':
            return ls_table[ sort_idx ]
        if self.sort_order == 'unsorted':
            return ls_table

    #---------------------------------
    # plot outputs in mlab-window 
    #---------------------------------

    plot_column = Enum( values = 'columns' )
    plot = Button
    def _plot_fired( self ):
        X = self.info_shell.X[:, 0]
        Y = self.info_shell.Y[:, 0]
        Z = self.info_shell.Z[:, 0]
        plot_col = getattr( self, self.plot_column )[:, 0]
        mlab.points3d( X, Y, Z, plot_col )
        mlab.show()

    #-------------------------------
    # ls group
    #-------------------------------

    # @todo: the dynamic selection of the columns to be displayed 
    # does not work in connection with the LSArrayAdapter 
    ls_group = VGroup( 
                        HGroup( Item( 'max_in_column' ),
                                Item( 'max_value', style = 'readonly', format_str = '%6.2f' ),
                                Item( 'max_value_all', style = 'readonly', format_str = '%6.2f' ),
                                Item( 'max_case', style = 'readonly', label = 'found in case: ' ),
                              ),
                        HGroup( Item( 'sort_column' ),
                                Item( 'sort_order' ),
                                Item( 'show_geo_columns', label = 'show geo' ),
                                Item( 'show_state_columns', label = 'show state' ),
                                Item( 'show_sr_columns', label = 'show sr' ),
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
    f_ctk = Float( 1.6, input = True )

    # flexural tensile strength [MPa]
    f_m = Float( 10.5, input = True )

    # ------------------------------------------------------------
    # SLS - derived params:
    # ------------------------------------------------------------

    # area
    #
    A = Property( Float, depends_on = 'info_shell.data_file_thickness' )
    def _get_A( self ):
        return self.info_shell.D_elem * 1.

    # moment of inertia
    #
    W = Property( Float, depends_on = 'info_shell.data_file_thickness' )
    def _get_W( self ):
        return 1. * self.info_shell.D_elem ** 2 / 6.

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

    #-------------------------------------------------------
    # get the maximum value and the corresponding case of 
    # the selected variable 'max_in_column' in all SLS sheets
    #-------------------------------------------------------

    max_value_all = Property( depends_on = 'max_in_column' )
    def _get_max_value_all( self ):
        return self.max_value_and_case['max_value']

    max_case = Property( depends_on = 'max_in_column' )
    def _get_max_case( self ):
        return self.max_value_and_case['max_case']

    max_value_and_case = Property( Dict, depends_on = 'max_in_column' )
    @cached_property
    def _get_max_value_and_case( self ):
        dir_list = DIRLIST
        sr_list = SRLIST
        ls_tree = self.info_shell.ls_tree
        max_value = 0.
        for dir in dir_list:
            for sr in sr_list:
                self.info_shell.ls_tree['SLS'][ sr ][ dir ].max_in_column = self.max_in_column
                max_value_ls = ls_tree['SLS'][ sr ][ dir ].max_value
                if max_value <= max_value_ls:
                    max_value = max_value_ls
                    max_case = 'S-' + sr + dir
        return { 'max_value':max_value,
                 'max_case':max_case }

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
                                Item( 'ls_table', show_label = False,
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

    # gamma-factor 
    gamma = Float( 1.5, input = True )

    # long term reduction factor
    beta = Float( 0.7, input = True )

    # INDEX l: longitudinal direction of the textile (MAG-02-02-06a)
    # characteristic tensile strength of the tensile specimen [N/mm2]
    f_tk_l = Float( 537, input = True )

    # design value of the tensile strength of the tensile specimen [N/mm2]
    # containing a gamma-factor of 1.5 and d long term reduction factor of 0.7
    # f_td_l = 251

    f_td_l = Property( Float, depends_on = '+input' )
    def _get_f_td_l( self ):
        return self.beta * self.f_tk_l / self.gamma

    # cross sectional area of the reinforcement [mm2/m]
    a_t_l = Float( 71.65, input = True )

    # INDEX q: orthogonal direction of the textile (MAG-02-02-06a)
    # characteristic tensile strength of the tensile specimen [N/mm2]
    f_tk_q = Float( 511, input = True )

    # design value of the tensile strength of the tensile specimen [kN/m]
    # f_td_q = 238
    f_td_q = Property( Float, depends_on = '+input' )
    def _get_f_td_q( self ):
        return self.beta * self.f_tk_q / self.gamma

    # cross sectional area of the reinforcement [mm2/m]
    a_t_q = Float( 53.31, input = True )

    # tensile strength of the textile reinforcement [kN/m]
    F_Rtex_l = Property( Float, depends_on = '+input' )
    def _get_F_Rtex_l( self ):
        return self.a_t_l * self.f_td_l / 1000.

    # tensile strength of the textile reinforcement [kN/m]
    F_Rtex_q = Property( Float )
    def _get_F_Rtex_q( self, depends_on = '+input' ):
        return self.a_t_q * self.f_td_q / 1000.

    # ------------------------------------------------------------
    # ULS - derived params:
    # ------------------------------------------------------------

    # Parameters for the cracked state (GdT):
    # assumptions!

    # (resultierende statische Nutzhoehe) 
    #
    d = Property( Float, depends_on = 'info_shell.data_file_thickness' )
    def _get_d( self ):
        return 0.75 * self.info_shell.D_elem

    # (Abstand Schwereachse zur resultierende Bewehrungslage) 
    # chose the same amount of reinforcement at the top as at the bottom 
    # i.e. zs = zs1 = zs2
    #
    zs = Property( Float, depends_on = 'info_shell.data_file_thickness' )
    def _get_zs( self ):
        return self.d - self.info_shell.D_elem / 2.

    # (Innerer Hebelarm) 
    #
    z = Property( Float )
    def _get_z( self ):
        return 0.9 * self.d

    # ------------------------------------------------------------
    # ULS: outputs
    # ------------------------------------------------------------

    ls_columns = List( [ 'e', 'm_Eds', 'f_t', 'beta_l', 'beta_q', 'f_Rtex', 'n_tex' ] )

#    sr_columns = [ 'm', 'n', 'alpha' ]
    sr_columns = [ 'm', 'n', 'alpha', 'd', 'zs', 'z' ]

    alpha_varname = Property()
    def _get_alpha_varname( self ):
        return 'alpha_' + self.stress_res

    alpha = Property
    def _get_alpha( self ):
        return getattr( self.info_shell, self.alpha_varname )

    ls_values = Property( depends_on = '+input' )
    @cached_property
    def _get_ls_values( self ):
        '''get the outputs for ULS
        '''
        n = self.n
        m = self.m
        alpha = self.alpha

        zs = self.zs
        z = self.z
        F_Rtex_l = self.F_Rtex_l
        F_Rtex_q = self.F_Rtex_q

        # (Exzentrizitaet)
        e = abs( m / n )
        e[ n == 0 ] = 1E9 # if normal force is zero set e to very large value

        # moment at the height of the resulting reinforcement layer:
        m_Eds = abs( m ) - zs * n

        # tensile force in the reinforcement for bending and compression
        f_t = m_Eds / z + n

        # check if the two conditions are true:
        cond1 = n > 0
        cond2 = e < zs
        bool_arr = cond1 * cond2
        # in case of pure tension in the cross section:
        f_t[ bool_arr ] = n[ bool_arr ] * ( zs[ bool_arr ] + e[ bool_arr ] ) / ( zs[ bool_arr ] + zs[ bool_arr ] )

        # angel of deflection of the textile reinforcement for dimensioning in x-direction
        # distinguished between longtudinal (l) and transversal (q) direction

        # ASSUMPTION: worst case angle used
        # as first step use the worst case reduction due to deflection possible (at 55 degrees)
        beta_l = 55. * pi / 180. * ones_like( alpha )
        beta_q = ( 90. - 55. ) * pi / 180. * ones_like( alpha )

        # @todo: as second step use the value for an alternating layup (i.e. deflection angle)
        # @todo: get the correct formula for the demonstrator arrangement 
        # i.e. the RFEM coordinate system orientation
#            beta_l = pi/2 - abs( alpha )
#            beta_q = abs( alpha )

        # resulting strength of the bi-directional textile considering the 
        # deflection of the reinforcement in the loading direction:
        f_Rtex = F_Rtex_l * cos( beta_l ) * ( 1 - beta_l / ( pi / 2 ) ) + \
                 F_Rtex_q * cos( beta_q ) * ( 1 - beta_q / ( pi / 2 ) )

        f_Rtex = 11.65 * ones_like( alpha )
        print 'NOTE: f_Rtex set to 11.65 kN/m !'

        # necessary number of reinfocement layers
        n_tex = f_t / f_Rtex

        return { 'e':e, 'm_Eds':m_Eds, 'f_t':f_t,
                 'beta_l':beta_l, 'beta_q':beta_q, 'f_Rtex':f_Rtex,
                 'n_tex':n_tex }

    e = Property
    def _get_e( self ):
        return self.ls_values['e']

    m_Eds = Property
    def _get_m_Eds( self ):
        return self.ls_values['m_Eds']

    f_t = Property
    def _get_f_t( self ):
        return self.ls_values['f_t']

    beta_l = Property
    def _get_beta_l( self ):
        return self.ls_values['beta_l']

    beta_q = Property
    def _get_beta_q( self ):
        return self.ls_values['beta_q']

    f_Rtex = Property
    def _get_f_Rtex( self ):
        return self.ls_values['f_Rtex']

    n_tex = Property
    def _get_n_tex( self ):
        return self.ls_values['n_tex']

    #-------------------------------------------------------
    # get the maximum value and the corresponding case of 
    # the selected variable 'max_in_column' in all ULS sheets
    #-------------------------------------------------------

    max_value_all = Property( depends_on = 'max_in_column' )
    def _get_max_value_all( self ):
        return self.max_value_and_case['max_value']

    max_case = Property( depends_on = 'max_in_column' )
    def _get_max_case( self ):
        return self.max_value_and_case['max_case']

    max_value_and_case = Property( Dict, depends_on = 'max_in_column' )
    @cached_property
    def _get_max_value_and_case( self ):
        dir_list = DIRLIST
        sr_list = SRLIST
        ls_tree = self.info_shell.ls_tree
        max_value = 0.
        for dir in dir_list:
            for sr in sr_list:
                self.info_shell.ls_tree['ULS'][ sr ][ dir ].max_in_column = self.max_in_column
                max_value_ls = ls_tree['ULS'][ sr ][ dir ].max_value
                if max_value <= max_value_ls:
                    max_value = max_value_ls
                    max_case = 'U-' + sr + dir
        return { 'max_value':max_value,
                 'max_case':max_case }

    #-------------------------------
    # ls view
    #-------------------------------

    # @todo: the dynamic selection of the columns to be displayed 
    # does not work in connection with the LSArrayAdapter 
    traits_view = View( 
                       VGroup( 
                        HGroup( 
                            VGroup( 
                                Item( name = 'gamma', label = 'security factor material [-]:  gamma ' ),
                                Item( name = 'beta', label = 'reduction long term durability [-]:  beta ' ),
                                label = 'security factors'
                                  ),
                            VGroup( 
                                Item( name = 'f_tk_l', label = 'characteristic strength textil [MPa]:  f_tk_l ', format_str = "%.1f" ),
                                Item( name = 'f_td_l', label = 'design strength textil [MPa]:  f_td_l ', style = 'readonly', format_str = "%.1f" ),
                                Item( name = 'a_t_l', label = 'cross sectional area textil [mm^2]:  a_t_l ', style = 'readonly', format_str = "%.1f" ),
                                Item( name = 'F_Rtex_l', label = 'Strength textil [kN/m]:  F_Rtex_l ', style = 'readonly', format_str = "%.0f" ),
                                label = 'material properties (longitudinal)'
                                  ),
                            VGroup( 
                                Item( name = 'f_tk_q', label = 'characteristic strength textil [MPa]:  f_tk_q ', format_str = "%.1f" ),
                                Item( name = 'f_td_q', label = 'design strength textil [MPa]:  f_td_q ', style = 'readonly', format_str = "%.1f" ),
                                Item( name = 'a_t_q', label = 'cross sectional area textil [mm^2]:  a_t_q ', style = 'readonly', format_str = "%.1f" ),
                                Item( name = 'F_Rtex_q', label = 'Strength textil [kN/m]: F_Rtex_q ', style = 'readonly', format_str = "%.0f" ),
                                label = 'material Properties (transversal)'
                                  ),
                             ),

                        VGroup( 
                            Include( 'ls_group' ),
                            Item( 'ls_table', show_label = False,
                                  editor = TabularEditor( adapter = LSArrayAdapter() ) )
                              ),
                            ),
                      resizable = True,
                      scrollable = True,
                      height = 1000,
                      width = 1100
                      )

LSLIST = [ SLS, ULS ]

class InfoShell( HasTraits ):
    '''Assessment tool
    '''

    lc_list = List( Instance( LC ), input = True )
    def _lc_list_default( self ):
        return [ LC( name = 'G', category = 'dead-load', file_name = 'input_data_stress_resultants.csv' ),

                 LC( name = 'G_A', category = 'additional dead-load', file_name = 'input_data_stress_resultants.csv' ),

                 LC( name = 'W (Druck)', category = 'imposed-load', file_name = 'input_data_stress_resultants.csv',
                     exclusive_to = ['W (Sog)'], psi_0 = 1.0, psi_1 = 1.0, psi_2 = 1.0 ),

                 LC( name = 'W (Sog)', category = 'imposed-load', file_name = 'input_data_stress_resultants.csv',
                     exclusive_to = ['W (Druck)'], psi_0 = 1.0, psi_1 = 1.0, psi_2 = 1.0 ),

                 LC( name = 'S', category = 'imposed-load', file_name = 'input_data_stress_resultants.csv',
                     exclusive_to = [], psi_0 = 1.0, psi_1 = 1.0, psi_2 = 1.0 )
               ]

    lcm = Property( Instance( LCManager ) )#, depends_on = 'lc_list' )
    def _get_lcm( self ):
        '''loading case manager'''
        return LCManager( lc_list = self.lc_list )

    lcc_tree = Property( Array )#, depends_on = 'lc_list' )
    def _get_lcc_tree( self ):
        '''loading case combination tree, e.g. lcc_tree['ULS']'''
        return self.lcm.lcc_tree

    current_ls = Enum( 'ULS', 'SLS' )

    #------------------------------------------
    # specify default data input files:
    #------------------------------------------

    # raw input file for thicknesses (and coordinates)
    #
    data_file_thickness = Str
    def _data_file_thickness_default( self ):
        return 'input_data_thickness.csv'

    #@todo: display all file_names and names of the lc's in the view

    # raw input file for element numbers, coordinates, and stress_resultants
    #
    data_file_stress_resultants = Str
    def _data_file_stress_resultants_default( self ):
        return self.lc_list[0].file_name

    #------------------------------------------
    # read the geometry data from file 
    # (corrds and thickness):
    #------------------------------------------

    def _read_thickness_data( self, file_name ):
        '''to read the stb-thickness save the xls-worksheet 
        to a csv-file using ';' as filed delimiter and ' ' (blank)
        as text delimiter.
        '''
        print '*** read thickness data ***'
        # get the column headings defined in the second row 
        # of the csv thickness input file
        # "Nr.;X;Y;Z;[mm]"
        #
        file = open( file_name, 'r' )
        first_line = file.readline()
        second_line = file.readline()
        column_headings = second_line.split( ';' )
        # remove '\n' from last string element in list
        column_headings[-1] = column_headings[-1][:-1]
        column_headings_arr = array( column_headings )
        elem_no_idx = where( 'Nr.' == column_headings_arr )[0]
        X_idx = where( 'X' == column_headings_arr )[0]
        Y_idx = where( 'Y' == column_headings_arr )[0]
        Z_idx = where( 'Z' == column_headings_arr )[0]
        thickness_idx = where( '[mm]' == column_headings_arr )[0]
        # read the float data:
        #
        input_arr = loadtxt( file_name, delimiter = ';', skiprows = 2 )

        # element number:
        #
        elem_no = input_arr[:, elem_no_idx]

        # coordinates [m]:
        #
        X_ = input_arr[:, X_idx]
        Y_ = input_arr[:, Y_idx]
        Z_ = input_arr[:, Z_idx]

        # element thickness [mm]:
        #
        thickness = input_arr[:, thickness_idx]

        return  {'X_':X_, 'Y_':Y_, 'Z_':Z_,
                 'thickness':thickness }

    # coordinates and element thickness read from file:
    # 
    thickness_data_dict = Property( Dict, depends_on = 'data_file_thickness' )
    @cached_property
    def _get_thickness_data_dict( self ):
        return self._read_thickness_data( self.data_file_thickness )

    X_ = Property( Array )
    def _get_X_( self ):
        return self.thickness_data_dict['X_']

    Y_ = Property( Array )
    def _get_Y_( self ):
        return self.thickness_data_dict['Y_']

    Z_ = Property( Array )
    def _get_Z_( self ):
        return self.thickness_data_dict['Z_']

    D_elem = Property( Array )
    def _get_D_elem( self ):
        '''element thickness (units changed form [mm] to [m])'''
        return self.thickness_data_dict['thickness'] / 1000.

    # ------------------------------------------------------------
    # Get the state data from the LC's 'state_data_dict' 
    # ------------------------------------------------------------

    # get elem_no and coordinates are taken from the first loading-case:
    # 
    state_data_dict = Property( Dict )#, depends_on = 'lc_list, current_ls' )
    @cached_property
    def _get_state_data_dict( self ):
        return self.lc_list[0].state_data_dict

    elem_no = Property( Array )
    def _get_elem_no( self ):
        return self.state_data_dict['elem_no']

    X = Property( Array )
    def _get_X( self ):
        return self.state_data_dict['X']

    Y = Property( Array )
    def _get_Y( self ):
        return self.state_data_dict['Y']

    Z = Property( Array )
    def _get_Z( self ):
        return self.state_data_dict['Z']

    # ------------------------------------------------------------
    # the stress resultants are taken from 'lcc_tree' depending 
    # on the current limit state  'current_ls' 
    # ------------------------------------------------------------

    mx = Property( Array, depends_on = 'current_ls' )
    @cached_property
    def _get_mx( self ):
        return self.lcc_tree[ self.current_ls ]['mx']

    my = Property( Array, depends_on = 'current_ls' )
    @cached_property
    def _get_my( self ):
        return self.lcc_tree[ self.current_ls ]['my']

    mxy = Property( Array, depends_on = 'current_ls' )
    @cached_property
    def _get_mxy( self ):
        return self.lcc_tree[ self.current_ls ]['mxy']

    nx = Property( Array, depends_on = 'current_ls' )
    @cached_property
    def _get_nx( self ):
        return self.lcc_tree[ self.current_ls ]['nx']

    ny = Property( Array, depends_on = 'current_ls' )
    @cached_property
    def _get_ny( self ):
        return self.lcc_tree[ self.current_ls ]['ny']

    nxy = Property( Array, depends_on = 'current_ls' )
    @cached_property
    def _get_nxy( self ):
        return self.lcc_tree[ self.current_ls ]['nxy']

    combi_key = Property( Array, depends_on = 'current_ls' )
    @cached_property
    def _get_combi_key( self ):
        return self.lcc_tree[ self.current_ls ]['combi_key']

    # ------------------------------------------------------------
    # check input files for consistency
    # ------------------------------------------------------------

    @on_trait_change( 'data_file_thickness, data_file_stress_resultants' )
    def _check_input_files_for_consistency( self ):
        '''Check if the coordinate order of the thickness input file is 
        identical to the order in the stress_resultant input file.

        Here the first sr-input file is taken. The internal consistency
        of all defined loading cases is checked in the LCManager.
        '''
        if not all( self.X ) == all( self.X_ ) or \
            not all( self.Y ) == all( self.Y_ ) or \
            not all( self.Z ) == all( self.Z_ ):
            raise ValueError, 'coordinates in file % s and file % s are not identical. Check input files for consistency ! ' \
                    % ( self.data_file_thickness, self.data_file_stress_resultants )
        else:
            print '*** input files checked for consistency ( OK ) *** '
            return True

    # ------------------------------------------------------------
    # Index M: calculate principle moments with corresponding normal forces
    # ------------------------------------------------------------

    princ_values_M = Property( Dict, depends_on = 'data_file_stress_resultants' )
    @cached_property
    def _get_princ_values_M( self ):
        '''principle value of the moments forces:
        and principle angle of the moments forces:
        mx_M, my_M, nx_M, ny_M: transform the values in the principle direction
        '''
        # stress_resultants in global coordinates
        #
        mx = self.mx
        my = self.my
        mxy = self.mxy
        nx = self.nx
        ny = self.ny
        nxy = self.nxy

        # principal values
        #
        m1 = 0.5 * ( mx + my ) + 0.5 * sqrt( ( mx - my ) ** 2 + 4 * mxy ** 2 )
        m2 = 0.5 * ( mx + my ) - 0.5 * sqrt( ( mx - my ) ** 2 + 4 * mxy ** 2 )
        alpha_M = pi / 2. * ones_like( m1 )
        bool = m2 != mx
        alpha_M[ bool ] = arctan( mxy[ bool ] / ( m2[ bool ] - mx[ bool ] ) )

        # transform to principal directions
        #
        mx_M = 0.5 * ( my + mx ) - 0.5 * ( my - mx ) * cos( 2 * alpha_M ) - mxy * sin( 2 * alpha_M )
        my_M = 0.5 * ( my + mx ) + 0.5 * ( my - mx ) * cos( 2 * alpha_M ) + mxy * sin( 2 * alpha_M )
        nx_M = 0.5 * ( ny + nx ) - 0.5 * ( ny - nx ) * cos( 2 * alpha_M ) - nxy * sin( 2 * alpha_M )
        ny_M = 0.5 * ( ny + nx ) + 0.5 * ( ny - nx ) * cos( 2 * alpha_M ) + nxy * sin( 2 * alpha_M )
        return { 'm1':m1, 'm2':m2, 'alpha_M':alpha_M,
                 'mx_M':mx_M, 'my_M':my_M,
                 'nx_M':nx_M, 'ny_M':ny_M }

    m1 = Property( Float )
    def _get_m1( self ):
        return self.princ_values_M['m1']

    m2 = Property( Float )
    def _get_m2( self ):
        return self.princ_values_M['m2']

    alpha_M = Property( Float )
    def _get_alpha_M( self ):
        return self.princ_values_M['alpha_M']

    mx_M = Property( Float )
    def _get_mx_M( self ):
        return self.princ_values_M['mx_M']

    my_M = Property( Float )
    def _get_my_M( self ):
        return self.princ_values_M['my_M']

    nx_M = Property( Float )
    def _get_nx_M( self ):
        return self.princ_values_M['nx_M']

    ny_M = Property( Float )
    def _get_ny_M( self ):
        return self.princ_values_M['ny_M']

    # ------------------------------------------------------------
    # Index N: principle normal forces with corresponding moments
    # ------------------------------------------------------------

    princ_values_N = Property( Dict, depends_on = 'data_file_stress_resultants' )
    @cached_property
    def _get_princ_values_N( self ):
        '''principle value of the normal forces:
        and principle angle of the normal forces:
        mx_N, my_N, nx_N, ny_N: transform the values in the principle normal direction
        '''
        # stress_resultants in global coordinates
        #
        mx = self.mx
        my = self.my
        mxy = self.mxy
        nx = self.nx
        ny = self.ny
        nxy = self.nxy

        # principal values
        #
        n1 = 0.5 * ( nx + ny ) + 0.5 * sqrt( ( nx - ny ) ** 2 + 4 * nxy ** 2 )
        n2 = 0.5 * ( nx + ny ) - 0.5 * sqrt( ( nx - ny ) ** 2 + 4 * nxy ** 2 )
        alpha_N = pi / 2. * ones_like( n1 )
        bool = n2 != nx
        alpha_N[ bool ] = arctan( nxy[ bool ] / ( n2[ bool ] - nx[ bool ] ) )
        # transform to principal directions
        mx_N = 0.5 * ( my + mx ) - 0.5 * ( my - mx ) * cos( 2 * alpha_N ) - mxy * sin( 2 * alpha_N )
        my_N = 0.5 * ( my + mx ) + 0.5 * ( my - mx ) * cos( 2 * alpha_N ) + mxy * sin( 2 * alpha_N )
        nx_N = 0.5 * ( ny + nx ) - 0.5 * ( ny - nx ) * cos( 2 * alpha_N ) - nxy * sin( 2 * alpha_N )
        ny_N = 0.5 * ( ny + nx ) + 0.5 * ( ny - nx ) * cos( 2 * alpha_N ) + nxy * sin( 2 * alpha_N )

        return{'n1' : n1, 'n2' : n2, 'alpha_N' : alpha_N,
               'mx_N' : mx_N, 'my_N' : my_N,
               'nx_N' : nx_N, 'ny_N' : ny_N }

    n1 = Property( Float )
    def _get_n1( self ):
        return self.princ_values_N['n1']

    n2 = Property( Float )
    def _get_n2( self ):
        return self.princ_values_N['n2']

    alpha_N = Property( Float )
    def _get_alpha_N( self ):
        return self.princ_values_N['alpha_N']

    mx_N = Property( Float )
    def _get_mx_N( self ):
        return self.princ_values_N['mx_N']

    my_N = Property( Float )
    def _get_my_N( self ):
        return self.princ_values_N['my_N']

    nx_N = Property( Float )
    def _get_nx_N( self ):
        return self.princ_values_N['nx_N']

    ny_N = Property( Float )
    def _get_ny_N( self ):
        return self.princ_values_N['ny_N']

    #------------------------------------------
    # combinations of limit states, stress resultants and directions
    #------------------------------------------

    ls_tree = Dict
    def _ls_tree_default( self ):
        dir_list = DIRLIST
        sr_list = SRLIST
        ls_list = LSLIST

        ls_dict = {}
        for ls_class in ls_list:
            sr_dict = {}
            for sr in sr_list:
                dir_dict = {}
                for dir in dir_list:
                    # set 'current_ls' to current 'ls_tree' ls-key
                    # the change of attribute 'current_ls' triggers 
                    # a 'depends_on' change and the Properties
                    # of the stress resultants are refreshed 
                    # 
                    current_ls = ls_class.__name__
                    self._reset_stress_resultants()
                    print 'current_ls', current_ls

                    dir_dict[ dir ] = ls_class( info_shell = self, dir = dir, sr = sr )
                sr_dict[ sr ] = dir_dict
            ls_dict[ ls_class.__name__ ] = sr_dict
        return ls_dict

#    for k, v in d.iteritems():
#        c[k] = v



    #------------------------------------------
    # get arrays for the TabularEditor:
    #------------------------------------------

    U_Mx = Property( Instance( ULS ) )
    def _get_U_Mx( self ):
        return self.ls_tree['ULS']['M']['x']

    U_My = Property( Instance( ULS ) )
    def _get_U_My( self ):
        return self.ls_tree['ULS']['M']['y']

    U_Nx = Property( Instance( ULS ) )
    def _get_U_Nx( self ):
        return self.ls_tree['ULS']['N']['x']

    U_Ny = Property( Instance( ULS ) )
    def _get_U_Ny( self ):
        return self.ls_tree['ULS']['N']['y']

    S_Mx = Property( Instance( SLS ) )
    def _get_S_Mx( self ):
        return self.ls_tree['SLS']['M']['x']

    S_My = Property( Instance( SLS ) )
    def _get_S_My( self ):
        return self.ls_tree['SLS']['M']['y']

    S_Nx = Property( Instance( SLS ) )
    def _get_S_Nx( self ):
        return self.ls_tree['SLS']['N']['x']

    S_Ny = Property( Instance( SLS ) )
    def _get_S_Ny( self ):
        return self.ls_tree['SLS']['N']['y']

    # ------------------------------------------------------------
    # View 
    # ------------------------------------------------------------

    traits_view = View( Item( 'data_file_stress_resultants', label = 'Evaluated input file for stress_resultants ',
                              style = 'readonly', emphasized = True ),

                        Item( 'data_file_thickness', label = 'Evaluated input file for thicknesses ',
                               style = 'readonly', emphasized = True ),

                        Item( 'S_Nx@' , label = "S-NX", show_label = False ),

#                        Tabbed( 
#                            Item( 'S_Nx@' , label = "S-NX", show_label = False ),
#                            Item( 'S_Ny@' , label = "S-NY", show_label = False ),
#                            Item( 'S_Mx@' , label = "S-MX", show_label = False ),
#                            Item( 'S_My@' , label = "S-MY", show_label = False ),
#                            Item( 'U_Nx@' , label = "U-NX", show_label = False ),
#                            Item( 'U_Ny@' , label = "U-NY", show_label = False ),
#                            Item( 'U_Mx@' , label = "U-MX", show_label = False ),
#                            Item( 'U_My@' , label = "U-MY", show_label = False ),
#                            scrollable = False,
#                         ),
                      resizable = True,
                      scrollable = True,
                      height = 1000,
                      width = 1100
                      )


if __name__ == '__main__':
    ifs = InfoShell()
    ifs.lcc_tree['ULS']['mx']
    ifs.configure_traits()


#    print ifs.columns
#
#    ifs.selected_dir = 'y'
#    print ifs.columns
#
#    ifs.selected_sr = 'N'
#    print ifs.columns
#
#    ifs.selected_ls = ULS
#    print ifs.columns
#
#    print 'n1'
#    print ifs.n1
#    print ifs.my_M
#
#    print 'n_table'
#    print ifs.current_ls.ls_table

