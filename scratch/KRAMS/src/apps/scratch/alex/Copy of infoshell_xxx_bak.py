
from enthought.traits.api import \
    HasTraits, Directory, List, Int, Float, Any, \
    on_trait_change, File, Constant, Instance, Trait, \
    Array, Str, Property, cached_property, WeakRef, \
    Dict, Button, Enum, Color

from enthought.util.home_directory import \
    get_home_directory

from enthought.traits.ui.api import \
    View, Item, DirectoryEditor, TabularEditor, HSplit, Tabbed, VGroup, \
    TableEditor, Group, ListEditor, VSplit, HSplit, HGroup, Spring


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



#-- global variables -------------------------------------------------

display_list_standard = ['elem_no',
                         'X',
                         'Y',
                         'D_elem',
                         'mx',
                         'my',
                         'mxy',
                         'nx',
                         'ny',
                         'nxy'
                         ]

display_list_withcaseidx = ['mx', 'my', 'nx', 'ny', 'alpha']

# order mus correspond to output order of 'get_outputs_GdG': 
#[ sig_n, sig_m, eta_n, eta_m, eta_tot ]
#
outputs_list_GdG = [ 'sig_n', 'sig_m', 'eta_n', 'eta_m', 'eta_tot' ]

# order mus correspond to output order of 'get_outputs_GdT': 
# [ e, m_Eds, f_t, beta_l, beta_q, f_Rtex, n_tex ]
#
outputs_list_GdT = ['e', 'm_Eds', 'f_t', 'beta_l', 'beta_q', 'f_Rtex', 'n_tex' ]

# -- case indices --------------
limit_state_list = [ 'G', 'T' ]
case_list = [ 'N', 'M' ]
direction_list = [ 'x', 'y' ]


#-- Tabular Adapter Definition -------------------------------------------------

class ArrayAdapter ( TabularAdapter ):

    columns = List()
    font = 'Courier 10'
    alignment = 'right'
    format = '%5.2f'#'%g'
    even_bg_color = Color( 0xE0E0FF )
    width = Float( 80 )

    # @temporary hack; is there a better solution?
    def get_format ( self, object, trait, row, column ):
        """ Returns the Python format string to use for a specified column.
        """
        if column in [0, 1, 2]:
            a = self.format
            b = '%3d'
            self.format = b
            c = self._result_for( 'get_format', object, trait, row, column )
            self.format = a
            return c
        else:
            return self._result_for( 'get_format', object, trait, row, column )

class ArrayAdapterGN ( ArrayAdapter ):
    '''print column headings in the first row
    columns = [ ( 'X', 0), ( 'Y', 1 ), ( 'Z', 2 ),
                ('mx', 3), ('my',4), ('mxy',5), ('nx',6), ('ny',7), ('nxy',8),\
                ('m1',9), ('m2',10), ('alpha_M',11),\
                ('mx_M',12), ('my_M',13), ('nx_M',14), ('ny_M',15),\
                ('sig_n',16), ('sig_m',17), ('eta_n',18), ('eta_m',19), ( 'eta_tot', 20 ) ]
    '''
    columns = []
    for idx, name in enumerate( display_list_standard ):
        columns.append( ( name, idx ) )
    for idx, name in enumerate( display_list_withcaseidx ):
        columns.append( ( name + '_N', idx ) )
    for idx, name in enumerate( outputs_list_GdG ):
        columns.append( ( name, idx ) )

class ArrayAdapterGM ( ArrayAdapter ):

    columns = []
    for idx, name in enumerate( display_list_standard ):
        columns.append( ( name, idx ) )

    for idx, name in enumerate( display_list_withcaseidx ):
        columns.append( ( name + '_M', idx ) )

    for idx, name in enumerate( outputs_list_GdG ):
        columns.append( ( name, idx ) )

class ArrayAdapterTN ( ArrayAdapter ):
    '''print column headings in the first row
    columns = [ ( 'X', 0), ( 'Y', 1 ), ( 'Z', 2 ), ('mx', 3), \
                ('my',4), ('mxy',5), ('nx',6), ('ny',7), ('nxy',8),\
                ('n1',9), ('n2',10), ('alpha_N',11),\
                ('mx_N',12), ('my_N',13), ('nx_N',14), ('ny_N',15),  \
                ('e',16), ('m_Eds',17), ('f_t',18), ( 'f_Rtex', 19 ), ( 'n_tex', 20 ) ]
    '''
    columns = []
    for idx, name in enumerate( display_list_standard ):
        columns.append( ( name, idx ) )
    for idx, name in enumerate( display_list_withcaseidx ):
        columns.append( ( name + '_N', idx ) )
    for idx, name in enumerate( outputs_list_GdT ):
        columns.append( ( name, idx ) )

class ArrayAdapterTM ( ArrayAdapter ):

    columns = []
    for idx, name in enumerate( display_list_standard ):
        columns.append( ( name, idx ) )
    for idx, name in enumerate( display_list_withcaseidx ):
        columns.append( ( name + '_M', idx ) )
    for idx, name in enumerate( outputs_list_GdT ):
        columns.append( ( name, idx ) )



class InfoShellFileMushroof( HasTraits ):
    '''
    Represent a single test specifying the design parameters.
    and access to the measured data.
    
    Made Assumptions for defined formulas:
    - reinforcement layers are spaced equally throughout the cross-section (d = 0,75*D_elem)
    - same reinforcement layers on top as on bottom of cross section (e = abs(m/n)
    '''

    # ------------------------------------------------------------
    # GdG: Input parameters for the uncracked state (G):
    # ------------------------------------------------------------

    # tensile strength [MPa]
    f_ctk = Float( 1.6, input = True )

    # flexural tensile strength [MPa]
    f_m = Float( 10.5, input = True )

    #--------------------------------------------------------
    # GdT: Input prameters for the cracked state / reinforcement 
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


    def __init__( self, **kw ):
        '''Initialization with traits.
        '''
        super( InfoShellFileMushroof, self ).__init__( **kw )

        for limit_state_name in limit_state_list:
            for case_name in case_list:
                for direction_name in direction_list:
                    self.add_trait( case_name + direction_name + '_' + limit_state_name, Array )

                    # get maximum loading case for the sort variable:
                    #
#                    for sort_var_name in [ self.sort_var_name_GdG, self.sort_var_name_GdT ]:
                    for sort_var_name in outputs_list_GdG + outputs_list_GdT:
                        self.add_trait( sort_var_name + '_max_case', Str )
                        self.add_trait( sort_var_name + '_max', Float )

                    for output_idx, output_name in enumerate( outputs_list_GdG + outputs_list_GdT ):
                        # assign: sig_n_Nx
                        #         sig_n_Ny
                        #         sig_n_Mx
                        #         sig_n_My
                        self.add_trait( output_name + '_' + case_name + direction_name, Array )

    #------------------------------------------
    # Data files and methods for file reading:
    #------------------------------------------

    # raw input file thickness
    #
    data_file_thickness = Str
    def _data_file_thickness_default( self ):
        return 'thickness_debug.csv'
        return 'thickness_mushroof_stb.csv'

    # raw input file solicitations
    #
    data_file_solicitations = Str
    def _data_file_solicitations_default( self ):
        return 'solicitations_debug.csv'
        return 'solicitations_mushroof_stb.csv'

    def read_data_file_thickness( self, file_name ):
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

        return  {'X_':X_,
                 'Y_':Y_,
                 'Z_':Z_,
                 'D_elem':thickness }

    def read_data_file_solicitations( self, file_name ):
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

        return { 'elem_no':elem_no,
                 'X':X,
                 'Y':Y,
                 'Z':Z,
                 'mx':mx,
                 'my':my,
                 'mxy':mxy,
                 'nx':nx,
                 'ny':ny,
                 'nxy':nxy }

    # ------------------------------------------------------------
    # Read input txt-file and assign input data 
    # ------------------------------------------------------------

    # coordinates and element thickness read from file:
    # 
    thickness_input_dict = Property( Dict, depends_on = 'data_file_thickness' )
    @cached_property
    def _get_thickness_input_dict( self ):
        return self.read_data_file_thickness( self.data_file_thickness )

    # coordinates and solicitations read from file:
    # 
    solicitations_input_dict = Property( Dict, depends_on = 'data_file_solicitations' )
    @cached_property
    def _get_solicitations_input_dict( self ):
        return self.read_data_file_solicitations( self.data_file_solicitations )

    # ------------------------------------------------------------
    # input attributes :
    # ------------------------------------------------------------

    def _assign_input_attributes( self ):
        '''assign the input data as array attributes to self
        '''
        # X_, Y_, Z_, D_elem 
        #
        for key, val in self.thickness_input_dict.iteritems():
            setattr( self, key, val )

        # elem_no, X, Y, Z, mx, my, mxy, nx, ny, nxy
        #
        for key, val in self.solicitations_input_dict.iteritems():
            setattr( self, key, val )

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

    princ_values_M_dict = Property( Array, depends_on = 'data_file_solicitations' )
    @cached_property
    def _get_princ_values_M_dict( self ):
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

    # ------------------------------------------------------------
    # Index N: principle normal forces with corresponding moments
    # ------------------------------------------------------------

    princ_values_N_dict = Property( Array, depends_on = 'data_file_solicitations' )
    @cached_property
    def _get_princ_values_N_dict( self ):
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

    # ------------------------------------------------------------
    # assign transformed solicitations:
    # ------------------------------------------------------------

    def _assign_derived_attributes( self ):
        '''assign the derived data as array attributes to self
        '''
        # e.g. mx_M, my_M, nx_M, ny_M
        #
        for key, val in self.princ_values_M_dict.iteritems():
            setattr( self, key, val )

        # mx_N, my_N, nx_N, ny_N
        #
        for key, val in self.princ_values_N_dict.iteritems():
            setattr( self, key, val )

    # ------------------------------------------------------------
    # GdG - derived params:
    # ------------------------------------------------------------

    # area
    #
    A = Property( Float, depends_on = 'data_file_thickness' )
    def _get_A( self ):
        return self.D_elem * 1.

    # moment of inertia
    #
    W = Property( Float, depends_on = 'data_file_thickness' )
    def _get_W( self ):
        return 1. * self.D_elem ** 2 / 6.

    # ------------------------------------------------------------
    # GdT - derived params:
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
    # GdG - outputs:
    # ------------------------------------------------------------ 
    # Index X: evaluation of stresses in transformed X-direction
    # Index Y: evaluation of stresses in transformed Y-direction

    def get_outputs_GdG_dict( self, n, m ):
        A = self.A
        W = self.W
        sig_n = n / A / 1000.
        sig_m = abs( m / W ) / 1000.
        eta_n = sig_n / self.f_ctk
        eta_m = sig_m / self.f_m
        eta_tot = eta_n + eta_m
        return { 'sig_n':sig_n,
                 'sig_m':sig_m,
                 'eta_n':eta_n,
                 'eta_m':eta_m,
                 'eta_tot':eta_tot }

    # ------------------------------------------------------------
    # T: outputs:
    # ------------------------------------------------------------

    def get_outputs_GdT_dict( self, n, m, alpha ):

        zs = self.zs
        z = self.z

        # (Exzentrizitaet)
        e = abs( m / n )
        e[n == 0] = 0 # if normal force is zero set e to zero

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
        f_Rtex = self.F_Rtex_l * cos( beta_l ) * ( 1 - beta_l / ( pi / 2 ) ) + \
                 self.F_Rtex_q * cos( beta_q ) * ( 1 - beta_q / ( pi / 2 ) )

        f_Rtex = 11.65 * ones_like( alpha )
        print 'NOTE: f_Rtex set to 11.65 kN/m !'

        # necessary number of reinfocement layers
        n_tex = f_t / f_Rtex

        return { 'e':e,
                 'm_Eds':m_Eds,
                 'f_t':f_t,
                 'beta_l':beta_l,
                 'beta_q':beta_q,
                 'f_Rtex':f_Rtex,
                 'n_tex':n_tex }


    # ------------------------------------------------------------
    # Output arrays for View
    # ------------------------------------------------------------

    # number of rows to be displayed
    result_arr_size = Int( 4 )

    ordering = Enum( "sorted", "unsorted" )

    #-------------------------------------------------------
    # variable to be used to determine the order of the results array
    #-------------------------------------------------------

    # e.g. 'eta_tot'
    sort_var_name_GdG = Enum( outputs_list_GdG[-1], outputs_list_GdG )
    # e.g. 'n_tex'
    sort_var_name_GdT = Enum( outputs_list_GdT[-1], outputs_list_GdT )

    #-------------------------------------------------------
    # assign outputs
    #-------------------------------------------------------

    def _assign_outputs( self ):
        '''assign outputs
        '''
        # e.g. value of 'eta_tot_max' for all cases '_Nx, _Ny, _Mx, _My'
        self.sort_var_max_list_GdG = []
        # e.g. value of 'n_tex_max' for all cases '_Nx, _Ny, _Mx, _My'
        self.sort_var_max_list_GdT = []

        for limit_state_name in limit_state_list:
            for case_name in case_list:
                for direction_name in direction_list:

                    #------------------------------------------
                    # solicitations
                    #------------------------------------------
                    # get: n = nx_N, m = mx_N
                    #      n = ny_N, m = my_N
                    #      n = nx_M, m = mx_M
                    #      n = ny_M, m = my_M
                    n = getattr( self, 'n' + direction_name + '_' + case_name )
                    m = getattr( self, 'm' + direction_name + '_' + case_name )
                    alpha = getattr( self, 'alpha' + '_' + case_name )

                    #------------------------------------------
                    # get outputs dict
                    #------------------------------------------

                    outputs_GdG_dict = self.get_outputs_GdG_dict( n, m )
                    outputs_GdT_dict = self.get_outputs_GdT_dict( n, m, alpha )

                    #------------------------------------------
                    # assign attributes
                    #------------------------------------------

                    # e.g. eta_tot_Nx¸ eta_tot_Ny, eta_tot_Mx, eta_tot_My
                    #
                    for key, val in outputs_GdG_dict.iteritems():
                        setattr( self, key, val )

                    # e.g. n_tex_Nx¸ n_tex_Ny, n_tex_Mx, n_tex_My
                    #
                    for key, val in outputs_GdT_dict.iteritems():
                        setattr( self, key, val )

    #-------------------------------------------------------
    # assign order index, collect max sort variables
    #-------------------------------------------------------

    def _assign_sort_index( self ):
        '''assign order index, collect max sort variables
        '''
        for limit_state_name in limit_state_list:
            for case_name in case_list:
                for direction_name in direction_list:

                    if limit_state_name == 'G':
                        sort_var_name = self.sort_var_name_GdG
                        print 'sort_var_name', sort_var_name

                    elif limit_state_name == 'T':
                        sort_var_name = self.sort_var_name_GdT
                        print 'sort_var_name', sort_var_name

                    else:
                        raise ValueError, 'limit_state_name does not exist'

                    # get eta_tot for the current case:
                    sort_var = getattr( self, sort_var_name + '_' + case_name + direction_name )

                    # get descending order of the values of 'eta_tot'

                    print 'SORTING ARRAY', sort_var

                    sort_var_idx = argsort( sort_var )[::-1]
                    print 'XXXX', argsort( sort_var )

                    print 'sort_var_idx', sort_var_idx[:10]
                    print 'sort_var_idx[0]', sort_var_idx[0]

                    # assign: Nx_G_idx 
                    #         Ny_G_idx 
                    #         Mx_G_idx 
                    #         My_G_idx 
                    setattr( self, case_name + direction_name + '_' + limit_state_name + '_idx', sort_var_idx )

                    # calculate the maximum value of eta_tot and n_tex 
                    # (displayed in the view as read-only parameters) 
                    #
                    sort_var_max = sort_var[ sort_var_idx[0] ]
#                    print 'sort_var_max', sort_var_max
#                    print 'sort_var_max (MAX)', max( sort_var )

                    # assign: eta_tot_Nx_max 
                    #         eta_tot_Ny_max 
                    #         eta_tot_Mx_max 
                    #         eta_tot_My_max 
                    setattr( self, sort_var_name + '_' + case_name + direction_name + '_max', sort_var_max )

                    # collect the case maximum in a list. This list is used to determin 
                    # the overall maximum of all cases in method '_assign_max_sort_var'
                    if limit_state_name == 'G':
#                        print 'append ', sort_var_max
                        self.sort_var_max_list_GdG.append( sort_var_max )
                    elif limit_state_name == 'T':
                        self.sort_var_max_list_GdT.append( sort_var_max )

                    # display values in the original order
                    #
                    if self.ordering == "unsorted":
                        print 'values unsorted'
                        idx_X = arange( self.elem_no.shape[0] )
                        setattr( self, case_name + direction_name + '_' + limit_state_name + '_idx', idx_X )


    #-------------------------------------------------------
    # Get max value of sort variable, number of displayed rows 
    #-------------------------------------------------------

    def _assign_max_sort_var( self ):
        '''Get max value of sort variable of all cases; 
        '''

        # get the max value of the sort variable and the corresponding loading case
        #
        self.sort_var_max_GdG = max( self.sort_var_max_list_GdG )
        self.sort_var_max_GdT = max( self.sort_var_max_list_GdT )

        for limit_state_name in limit_state_list:
            for case_name in case_list:
                for direction_name in direction_list:

                    # get maximum loading case for the sort variable:
                    #
                    if limit_state_name == 'G':
                        sort_var_name = self.sort_var_name_GdG
                        sort_var_max = getattr( self, sort_var_name + '_' + case_name + direction_name + '_max' )
                        if sort_var_max == self.sort_var_max_GdG:
#                            print 'max case found GdG'
#                            print 'sort_var_name', sort_var_name
#                            print 'sort_var_max', sort_var_max
#                            print 'sort_var_max_GdG', sort_var_max_GdG
                            setattr( self, sort_var_name + '_max_case', 'G-' + case_name + direction_name )
                            setattr( self, sort_var_name + '_max', sort_var_max )
                    elif limit_state_name == 'T':
                        sort_var_name = self.sort_var_name_GdT
                        sort_var_max = getattr( self, sort_var_name + '_' + case_name + direction_name + '_max' )
#                        print 'sort_var_name', sort_var_name
#                        print 'sort_var_max', sort_var_max
#                        print 'sort_var_max_GdT', self.sort_var_max_GdT

                        if sort_var_max == self.sort_var_max_GdT:
                            print 'max case found GdT'
                            setattr( self, sort_var_name + '_max_case', 'T-' + case_name + direction_name )
                            setattr( self, sort_var_name + '_max', sort_var_max )


    display_array = Array( float )
    # ------------------------------------------------------------
    # stack attributes together to form output arrays used for table View 
    # ------------------------------------------------------------
    def _assign_display_arrays( self ):
        '''stack attributes together for output arrays used for table View
        '''
        # number of displayed results
        #
        NRES = self.result_arr_size

        # get the displayed array for attributes defined in 'display_list'
        #
        for limit_state_name in limit_state_list:
            for case_name in case_list:
                for direction_name in direction_list:

                    # first columns (without case index)
                    display_arr = getattr( self, display_list_standard[0] )[:, newaxis]
                    for display_name in display_list_standard[1:]:
                        col = getattr( self, display_name )[:, newaxis]
                        display_arr = hstack( [ display_arr, col ] )

                    # additional column (with case index, e.g. '_Nx')
                    for display_name in display_list_withcaseidx:
                        col = getattr( self, display_name + '_' + case_name )[:, newaxis]
                        display_arr = hstack( [ display_arr, col ] )

                    # additional column (without case index)
                    if limit_state_name == 'G':
                        outputs_list = outputs_list_GdG
                    elif limit_state_name == 'T':
                        outputs_list = outputs_list_GdT
                    for display_name in outputs_list:
                        col = getattr( self, display_name + '_' + case_name + direction_name )[:, newaxis]
                        display_arr = hstack( [ display_arr, col ] )

                    # sort the values based on the order index
                    idx_X = getattr( self, case_name + direction_name + '_' + limit_state_name + '_idx' )
#                    print 'idx.shape', idx_X.shape
#                    print 'idx', self.Ny_T_idx

                    # assign: Nx_G, Ny_G, Mx_G, My_G
                    #
                    trait_name = case_name + direction_name + '_' + limit_state_name
                    arr = display_arr[ idx_X, : ][:NRES, :]
                    setattr( self, trait_name, arr )

    # ------------------------------------------------------------
    # define column headings
    # ------------------------------------------------------------

    columns_GN = []
    for idx, name in enumerate( display_list_standard ):
        columns_GN.append( ( name, idx ) )
    for idx, name in enumerate( display_list_withcaseidx ):
        columns_GN.append( ( name + '_N', idx ) )
    for idx, name in enumerate( outputs_list_GdG ):
        columns_GN.append( ( name, idx ) )

    # ------------------------------------------------------------
    # run evaluation
    # ------------------------------------------------------------

    @on_trait_change( '+input, data_file_thickness, data_file_solicitations, result_arr_size, f_ctk, f_m, ordering,\
                        gamma, beta, f_tk_l, f_tk_q, sort_var_name_GdG, sort_var_name_GdT' )
    def _run_evaluation( self ):

        print '*** run evaluation ***'

        self._assign_input_attributes()
        self._assign_derived_attributes()
        self._check_input_files_for_consistency()
        self._assign_outputs()
        self._assign_sort_index()
        self._assign_max_sort_var()
        self._assign_display_arrays()
#        print 'e', self.e_Nx
#        print 'n_tex', self.n_tex_Nx
#        print 'display_arr8', self.Nx_G
#        print 'self.sort_var_max_list_GdT', self.sort_var_max_list_GdT
#        print 'self.sort_var_max_GdT', self.sort_var_max_GdT
#        print 'n_tex_Nx (unordered)', self.n_tex_Nx
#        print 'n_tex_Nx (ordered)', self.n_tex_Nx[ self.Nx_T_idx]
#        print 'Nx_T[-1]', self.Nx_T[:, -1][ argsort( self.Nx_T[:, -1][::-1] )]
#        print 'Nx_T[-2]', self.Nx_T[:, -2]
#        print 'f_Rtex_N', self.f_Rtex_Nx
#        print 'display_arr', self.display_arr

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
                        HSplit( 
                        VGroup( Item( name = 'sort_var_name_GdG', label = 'sort_var_name' ),
                               Item( name = 'eta_tot_max', label = "Maximum value ", style = 'readonly' ),
                               Item( name = 'eta_tot_max_case', label = 'reached in case', style = 'readonly' ),
                               label = 'GdG',
                               dock = 'tab'
                               ),
                        VGroup( Item( name = 'sort_var_name_GdT', label = 'sort_var_name' ),
                               Item( name = 'n_tex_max', label = "Maximum value ", style = 'readonly' ),
                               Item( name = 'n_tex_max_case', label = 'reached in case', style = 'readonly' ),
                               label = 'GdT',
                               dock = 'tab'
                               ),
                        VGroup( Item( 'ordering', label = 'Use ordering', style = 'custom' ),
                               Item( name = 'result_arr_size', label = 'Number of selected rows' ),
                               label = 'display_options',
                               dock = 'tab'
                               ),
                               ),
#                        HSplit(Item( name = sort_var_name_GdG + '_max', label='Maximum value of ' + sort_var_name_GdG + '_max        ', style = 'readonly' ),
#                               Item( name = sort_var_name_GdG + '_max_case', label = 'reached in case', style = 'readonly' )),
#                        HSplit(Item( name = sort_var_name_GdT + '_max', label='Maximum value of ' + sort_var_name_GdT + '_max        ', style = 'readonly' ),
#                               Item( name = sort_var_name_GdT + '_max_case', label = 'reached in case', style = 'readonly' )),

                        Tabbed( 
                            Item( 'Nx_G' , label = "G-NX", show_label = False, editor = TabularEditor( adapter = ArrayAdapter( columns = columns_GN ) ) ),
                            Item( 'Ny_G' , label = "G-NY", show_label = False, editor = TabularEditor( adapter = ArrayAdapterGN() ) ),
                            Item( 'Mx_G' , label = "G-MX", show_label = False, editor = TabularEditor( adapter = ArrayAdapterGM() ) ),
                            Item( 'My_G' , label = "G-MY", show_label = False, editor = TabularEditor( adapter = ArrayAdapterGM() ) ),
                            Item( 'Nx_T' , label = "T-NX", show_label = False, editor = TabularEditor( adapter = ArrayAdapterTN() ) ),
                            Item( 'Ny_T' , label = "T-NY", show_label = False, editor = TabularEditor( adapter = ArrayAdapterTN() ) ),
                            Item( 'Mx_T' , label = "T-MX", show_label = False, editor = TabularEditor( adapter = ArrayAdapterTM() ) ),
                            Item( 'My_T' , label = "T-MY", show_label = False, editor = TabularEditor( adapter = ArrayAdapterTM() ) ),
                            scrollable = False,
                         ),
                      resizable = True,
                      scrollable = True,
                      height = 1000,
                      width = 1100
                      )


if __name__ == '__main__':
    isf = InfoShellFileMushroof( data_file_thickness = 'thickness_debug.csv',
                                 data_file_solicitations = 'solicitations_debug.csv' )
#    isf = InfoShellFileMushroof( data_file_thickness = 'thickness_mushroof_stb.csv',
#                                 data_file_solicitations = 'solicitations_mushroof_stb.csv' )
    isf._run_evaluation()
    isf.configure_traits()
