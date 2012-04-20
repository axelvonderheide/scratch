'''
Created on Jun 23, 2010

@author: alexander
'''

from enthought.traits.api import \
    HasTraits, Directory, List, Int, Float, Any, Enum, \
    on_trait_change, File, Constant, Instance, Trait, \
    Array, Str, Property, cached_property, WeakRef, \
    Dict, Button, Color

from enthought.util.home_directory import \
    get_home_directory

from enthought.traits.ui.api import \
    View, Item, DirectoryEditor, TabularEditor, HSplit, Tabbed, VGroup, \
    TableEditor, Group, ListEditor, VSplit, HSplit, HGroup, Spring

from enthought.mayavi import mlab

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

#class ArrayAdapter ( TabularAdapter ):
#    '''Array adapter for TabularEditor '''
#    current_ls = Instance( InfoShell )
#
#    columns = Property( List )
#    def _get_columns( self ):
#        return self.current_ls.columns



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

#    adapter_column_map = Property(depends_on = 'adapters,columns')



#    # @temporary hack; is there a better solution?
#    def get_format ( self, object, trait, row, column ):
#        """ Returns the Python format string to use for a specified column.
#        """
#        if column in [0, 1, 2]:
#            a = self.format
#            b = '%3d'
#            self.format = b
#            c = self._result_for( 'get_format', object, trait, row, column )
#            self.format = a
#            return c
#        else:
#            return self._result_for( 'get_format', object, trait, row, column )

class LS( HasTraits ):
    '''Limit state class
    '''

    # backward link to the info shell to access the
    # input data when calculating 
    # the limit-state-specific values
    #
    info_shell = WeakRef

    # parameters of the limit state

    dir = Enum( DIRLIST )

    stress_res = Enum( SRLIST )

    # limit state columns - defined in the subclasses
    ls_columns = List

    # stress resultant columns - generated from the parameter combination
    # dir and stress_res - one of MX, NX, MY, NY
    # e.g. [ mx_N, nx_N ] of [ my_M, ny_M], etc. 
    #
    m_varname = Property( Str )
    def _get_m_varname( self ):
        appendix = self.dir + '_' + self.stress_res
        return 'm' + appendix

    n_varname = Property( Str )
    def _get_n_varname( self ):
        appendix = self.dir + '_' + self.stress_res
        return 'n' + appendix

    n = Property( Float )
    def _get_n( self ):
        return getattr( self.info_shell, self.n_varname )

    m = Property( Float )
    def _get_m( self ):
        return getattr( self.info_shell, self.m_varname )

    sr_columns = List( ['m', 'n'] )

    #-------------------------------
    # geo columns form info shell
    #-------------------------------

    geo_columns = Property( Int )
    def _get_geo_columns( self ):
        return self.info_shell.geo_columns

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

    #-------------------------------
    # state columns form info shell
    #-------------------------------

    geo_columns = List( [ 'elem_no', 'X', 'Y', 'Z' ] )

    state_columns = List( ['mx', 'my', 'mxy', 'nx', 'ny', 'nxy' ] )


    state_columns = Property( Int )
    def _get_state_columns( self ):
        return self.info_shell.state_columns

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

    # all columns associated with the limit state including the corresponding
    # stress resultants
    #
    columns = Property
    def _get_columns( self ):
        return self.geo_columns + self.sr_columns + self.ls_columns

    sort_column = Enum( values = 'columns' )

    ls_table = Property( Array, depends_on = 'sort_column' )
    @cached_property
    def _get_ls_table( self ):

        arr_list = [ getattr( self, col ) for col in self.columns ]

        # get the array currently selected by the sort_column enumeration
        #
        sort_arr = getattr( self, self.sort_column )[:, 0]
        sort_idx = argsort( sort_arr )
        ls_table = hstack( arr_list )
        return ls_table[ sort_idx[::-1] ]


    plot_column = Enum( values = 'columns' )
    plot = Button
    def _plot_fired( self ):
        X = self.info_shell.X[:, 0]
        Y = self.info_shell.Y[:, 0]
        Z = self.info_shell.Z[:, 0]
        plot_col = getattr( self, self.plot_column )[:, 0]
        mlab.points3d( X, Y, Z, plot_col )
        mlab.show()

    def default_traits_view( self ):
        return View( HGroup( Item( 'sort_column' ),
                             Item( 'plot_column' ),
                             Item( 'plot' ), ),
                     Item( 'ls_table' ,
                            show_label = False,
                            editor = TabularEditor( adapter = LSArrayAdapter() ) ),
                      resizable = True,
                      scrollable = True,
                      height = 1000,
                      width = 1100
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
    A = Property( Float, depends_on = 'infoshell.data_file_thickness' )
    def _get_A( self ):
        return self.info_shell.D_elem * 1.

    # moment of inertia
    #
    W = Property( Float, depends_on = 'infoshell.data_file_thickness' )
    def _get_W( self ):
        return 1. * self.info_shell.D_elem ** 2 / 6.

    # ------------------------------------------------------------
    # SLS: outputs
    # ------------------------------------------------------------

    ls_columns = List( ['sig_n', 'sig_m', 'eta_n', 'eta_m', 'eta_tot', ] )

    ls_values = Property( depends_on = '+input' )
    @cached_property
    def _get_ls_values( self ):

        n = self.n
        m = self.m
        A = self.A
        W = self.W
        f_ctk = self.info_shell.f_ctk
        f_m = self.info_shell.f_m

        sig_n = n / A / 1000.
        sig_m = abs( m / W ) / 1000.
        eta_n = sig_n / f_ctk
        eta_m = sig_m / f_m
        eta_tot = eta_n + eta_m

        return { 'sig_n':sig_n,
                 'sig_m':sig_m,
                 'eta_n':eta_n,
                 'eta_m':eta_m,
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
    d = Property( Float, depends_on = 'data_file_thickness' )
    def _get_d( self ):
        return 0.75 * self.D_elem

    # (Abstand Schwereachse zur resultierende Bewehrungslage) 
    # chose the same amount of reinforcement at the top as at the bottom 
    # i.e. zs = zs1 = zs2
    #
    zs = Property( Float, depends_on = 'data_file_thickness' )
    def _get_zs( self ):
        return self.d - self.D_elem / 2.

    # (Innerer Hebelarm) 
    #
    z = Property( Float, depends_on = 'data_file_thickness' )
    def _get_z( self ):
        return 0.9 * self.d

    # ------------------------------------------------------------
    # ULS: outputs
    # ------------------------------------------------------------

    ls_columns = List( [ 'e', 'm_Eds', 'f_t', 'beta_l', 'beta_q', 'f_Rtex', 'n_tex' ] )

    sr_columns = ['m', 'n', 'alpha']

    alpha_varname = Property()
    def _get_alpha_varname( self ):
        return 'alpha_' + self.stress_res

    alpha = Property
    def _get_alpha( self ):
        return getattr( self.info_shell, self.alpha_varname )

    ls_values = Property( depends_on = 'info_shell.input' )
    @cached_property
    def _get_ls_values( self ):

        n = self.n
        m = self.m
        alpha = self.alpha
        zs = self.info_shell.zs
        z = self.info_shell.z
        # (Exzentrizitaet)
        e = abs( m / n )
        e[ n == 0 ] = 1E9 # if normal force is zero set e to very large value
        # @todo: verify formulas

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
        f_Rtex = self.info_shell.F_Rtex_l * cos( beta_l ) * ( 1 - beta_l / ( pi / 2 ) ) + \
                 self.info_shell.F_Rtex_q * cos( beta_q ) * ( 1 - beta_q / ( pi / 2 ) )

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


LSLIST = [ SLS, ULS ]

class InfoShell( HasTraits ):
    '''Assessment tool
    '''
    #------------------------------------------
    # specify default data files:
    #------------------------------------------

    # raw input file thickness
    #
    data_file_thickness = Str
    def _data_file_thickness_default( self ):
#        return 'thickness_debug.csv'
        return 'thickness_mushroof_stb.csv'

    # raw input file solicitations
    #
    data_file_solicitations = Str
    def _data_file_solicitations_default( self ):
#        return 'solicitations_debug.csv'
        return 'solicitations_mushroof_stb.csv'

    #------------------------------------------
    # read the geometry data from file 
    # (corrds and thickness):
    #------------------------------------------

    def _read_thickness_data( self, file_name ):
        '''to read the stb-thickness save the xls-worksheet 
        to a csv-file using ';' as filed delimiter and ' ' (blank)
        as text delimiter.
        '''
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

        # element thickness (units changed form [mm] to [m]):
        #
        thickness = input_arr[:, thickness_idx] / 1000.

        return  {'X_':X_, 'Y_':Y_, 'Z_':Z_,
                 'D_elem':thickness }

    # coordinates and element thickness read from file:
    # 
    thickness_data = Property( Dict, depends_on = 'data_file_thickness' )
    @cached_property
    def _get_thickness_data( self ):
        return self._read_thickness_data( self.data_file_thickness )

    X_ = Property( Array )
    def _get_X_( self ):
        return self.thickness_data['X_']

    Y_ = Property( Array )
    def _get_Y_( self ):
        return self.thickness_data['Y_']

    Z_ = Property( Array )
    def _get_Z_( self ):
        return self.thickness_data['Z_']

    D_elem = Property( Array )
    def _get_D_elem( self ):
        # change units from [mm] to [m]:
        return self.thickness_data['D_elem'] / 1000.0

    #------------------------------------------
    # read the state data from file 
    # (elem_no, coords, moments and normal forces):
    #------------------------------------------

    def _read_state_data( self, file_name ):
        '''to read the stb-solicitations save the xls-worksheet 
        to a csv-file using ';' as filed delimiter and ' ' (blank)
        as text delimiter.
        '''

        # this method returns only the MAXIMUM VALUES!!!
        # @todo: dublicate the elem number and coordinates and add also the minimum values

        # get the column headings defined in the second row 
        # of the csv soliciotations input file
        #
#        column_headings = array(["Nr.","Punkt","X","Y","Z","mx","my","mxy","vx","vy","nx","ny","nxy"])
        file = open( file_name, 'r' )
        lines = file.readlines()
        column_headings = lines[1].split( ';' )
        # remove '\n' from last string element in list
        column_headings[-1] = column_headings[-1][:-1]
        column_headings_arr = array( column_headings )
        elem_no_idx = where( 'Nr.' == column_headings_arr )[0]
        X_idx = where( 'X' == column_headings_arr )[0]
        Y_idx = where( 'Y' == column_headings_arr )[0]
        Z_idx = where( 'Z' == column_headings_arr )[0]
        mx_idx = where( 'mx' == column_headings_arr )[0]
        my_idx = where( 'my' == column_headings_arr )[0]
        mxy_idx = where( 'mxy' == column_headings_arr )[0]
        nx_idx = where( 'nx' == column_headings_arr )[0]
        ny_idx = where( 'ny' == column_headings_arr )[0]
        nxy_idx = where( 'nxy' == column_headings_arr )[0]

        # define arrays containing the information from the raw input file
        #

        # @todo: check how loadtxt can be used directly
        #        instead of reading lines line by line?
        # input_arr = loadtxt( file_name , delimiter=';', skiprows = 2 )

        # read max-values (=first row of each double-line):
        # read solicitation csv-input file line by line in steps of two
        # starting in the third line returning an data array.
        #
        n_columns = len( lines[2].split( ';' ) )
        # auxiliary line in array sliced off in input array
        data_array = zeros( n_columns )
        for n in range( 2, len( lines ), 2 ):
            line_split = lines[n].split( ';' )
            line_array = array( line_split, dtype = float )
            data_array = vstack( [ data_array, line_array ] )
        input_arr = data_array[1:, :]

        # element number:
        #
        elem_no = input_arr[:, elem_no_idx]

        # coordinates [m]:
        #
        X = input_arr[:, X_idx]
        Y = input_arr[:, Y_idx]
        Z = input_arr[:, Z_idx]
        coords_xyz = c_[ X, Y, Z ]

        # moments [kNm/m]
        #
        mx = input_arr[:, mx_idx]
        my = input_arr[:, my_idx]
        mxy = input_arr[:, mxy_idx]

        # normal forces [kN/m]:
        #
        nx = input_arr[:, nx_idx]
        ny = input_arr[:, ny_idx]
        nxy = input_arr[:, nxy_idx]

        return { 'elem_no':elem_no, 'X':X, 'Y':Y, 'Z':Z,
                 'mx':mx, 'my':my, 'mxy':mxy,
                 'nx':nx, 'ny':ny, 'nxy':nxy }

    # ------------------------------------------------------------
    # Read state data from file and assign attributes 
    # ------------------------------------------------------------

    # coordinates and solicitations read from file:
    # 
    state_data = Property( Dict, depends_on = 'data_file_solicitations' )
    @cached_property
    def _get_state_data( self ):
        return self._read_state_data( self.data_file_solicitations )

    elem_no = Property( Array )
    def _get_elem_no( self ):
        return self.state_data['elem_no']

    X = Property( Array )
    def _get_X( self ):
        return self.state_data['X']

    Y = Property( Array )
    def _get_Y( self ):
        return self.state_data['Y']

    Z = Property( Array )
    def _get_Z( self ):
        return self.state_data['Z']

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
    # check input files for consistency
    # ------------------------------------------------------------

    def _check_input_files_for_consistency( self ):
        '''check if the element order of the thickness input file is 
        identical to the order in the solicitation input file
        '''
        if any( 1 - ( self.X / self.X_ ) ) > 0.02 and \
            any( 1 - ( self.Y / self.Y_ ) ) > 0.02 and \
            any( 1 - ( self.Z / self.Z_ ) ) > 0.02:
            raise ValueError, 'coordinates in file %s and file %s are not identical. Check input files for consistency!' \
                    % ( self.data_file_thickness, self.data_file_solicitations )
        else:
            print '*** input files checked for consistency (OK) ***'


    # ------------------------------------------------------------
    # Index M: calculate principle moments with corresponding normal forces
    # ------------------------------------------------------------

    princ_values_M = Property( Dict, depends_on = 'data_file_solicitations' )
    @cached_property
    def _get_princ_values_M( self ):
        '''principle value of the moments forces:
        and principle angle of the moments forces:
        mx_M, my_M, nx_M, ny_M: transform the values in the principle direction
        '''
        # solicitations in global coordinates
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
        return { 'm1':m1,
                 'm2':m2,
                 'alpha_M':alpha_M,
                 'mx_M':mx_M,
                 'my_M':my_M,
                 'nx_M':nx_M,
                 'ny_M':ny_M }

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

    princ_values_N = Property( Dict, depends_on = 'data_file_solicitations' )
    @cached_property
    def _get_princ_values_N( self ):
        '''principle value of the normal forces:
        and principle angle of the normal forces:
        mx_N, my_N, nx_N, ny_N: transform the values in the principle normal direction
        '''
        # solicitations in global coordinates
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

        return{'n1' : n1,
               'n2' : n2,
               'alpha_N' : alpha_N,
               'mx_N' : mx_N,
               'my_N' : my_N,
               'nx_N' : nx_N,
               'ny_N' : ny_N }

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

    selected_ls = Enum( LSLIST )
    selected_dir = Enum( DIRLIST )
    selected_sr = Enum( SRLIST )

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
                    dir_dict[ dir ] = ls_class( info_shell = self, dir = dir, sr = sr )
                sr_dict[ sr ] = dir_dict
            ls_dict[ ls_class.__name__ ] = sr_dict
        return ls_dict

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

    traits_view = View( Item( 'data_file_solicitations', label = 'Evaluated input file for solicitations ',
                              style = 'readonly', emphasized = True ),

                        Item( 'data_file_thickness', label = 'Evaluated input file for thicknesses ',
                               style = 'readonly', emphasized = True ),

                        Item( name = 'f_ctk', label = 'Tensile strength concrete [MPa]:  f_ctk      ' ),
                        Item( name = 'f_m', label = 'Flexural tensile trength concrete [MPa]:  f_m        ' ),

                        Item( name = 'gamma', label = 'Security factor material [-]:  gamma ' ),
                        Item( name = 'beta', label = 'reduction factor for long term durability [-]:  beta ' ),

                        HSplit( 
                            Spring(),
                            VGroup( 
                                Item( name = 'f_tk_l', label = 'characteristic strength textil (longitudinal) [MPa]:  f_tk_l ', format_str = "%.1f" ),
                                Item( name = 'f_td_l', label = 'design strength textil (longitudinal) [MPa]:  f_td_l ', style = 'readonly', format_str = "%.1f" ),
                                Item( name = 'a_t_l', label = 'cross sectional area textil (longitudinal) [mm^2]:  a_t_l ', style = 'readonly', format_str = "%.1f" ),
                                Item( name = 'F_Rtex_l', label = 'Strength textil (longitudinal) [kN/m]:  F_Rtex_l ', style = 'readonly', format_str = "%.0f" ),
                                  ),
                            Spring(),
                            VGroup( 
                                Item( name = 'f_tk_q', label = 'characteristic strength textil (transversal) [MPa]:  f_tk_q ', format_str = "%.1f" ),
                                Item( name = 'f_td_q', label = 'design strength textil (transversal) [MPa]:  f_td_q ', style = 'readonly', format_str = "%.1f" ),
                                Item( name = 'a_t_q', label = 'cross sectional area textil (transversal) [mm^2]:  a_t_q ', style = 'readonly', format_str = "%.1f" ),
                                Item( name = 'F_Rtex_q', label = 'Strength textil (transversal) [kN/m]: F_Rtex_q ', style = 'readonly', format_str = "%.0f" ),
                                  ),
                            Spring(),
                             ),
                        Tabbed( 
                            Item( 'S_Nx@' , label = "G-NX", show_label = False ),
                            Item( 'S_Ny@' , label = "G-NY", show_label = False ),
                            Item( 'S_Mx@' , label = "G-MX", show_label = False ),
                            Item( 'S_My@' , label = "G-MY", show_label = False ),
                            Item( 'U_Nx@' , label = "T-NX", show_label = False ),
                            Item( 'U_Ny@' , label = "T-NY", show_label = False ),
                            Item( 'U_Mx@' , label = "T-MX", show_label = False ),
                            Item( 'U_My@' , label = "T-MY", show_label = False ),
                            scrollable = False,
                         ),
                      resizable = True,
                      scrollable = True,
                      height = 1000,
                      width = 1100
                      )



if __name__ == '__main__':
    ifs = InfoShell()

    #print ifs.columns

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

    ifs.configure_traits()

