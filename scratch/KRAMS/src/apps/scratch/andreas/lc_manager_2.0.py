'''
Created on Jun 29, 2010

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

from numpy import \
    array, loadtxt, arange, sqrt, zeros, arctan, sin, cos, ones_like, \
    vstack, savetxt, hstack, argsort, fromstring, zeros_like, \
    copy, c_, newaxis, argmax, where, sqrt, frompyfunc, sum, \
    ones, transpose, shape, append

from itertools import \
    product

from math import pi
from string import split
import os

from scipy.io import read_array

class LC( HasTraits ):
    '''Loading case class
    '''
    # name of the file containing the stress resultants
    # of the loading case
    #
    file_name = Str

    # name of the loading case
    #
    name = Str

    # category of the loading case
    #
    category = Enum( 'dead-load', 'additional dead-load', 'imposed-load' )

    # list of keys specifing other loading cases 
    # that can not exist at the same time
    # 
    exclusive_to = List( Str )
    def _exclusive_to_default( self ):
        return []

    # security factors
    #
    gamma_fav = Float
    def _gamma_fav_default( self ):
        if self.category == 'dead-load':
            return 1.00
        if self.category == 'additional dead-load':
            return 0.00
        if self.category == 'imposed-load':
            return 0.00

    gamma_unf = Float
    def _gamma_unf_default( self ):
        if self.category == 'dead-load':
            return 1.35
        if self.category == 'additional dead-load':
            return 1.35
        if self.category == 'imposed-load':
            return 1.50

    gamma_fav_SLS = Float
    def _gamma_fav_SLS_default( self ):
        if self.category == 'dead-load':
            return 1.00
        if self.category == 'additional dead-load':
            return 0.00
        if self.category == 'imposed-load':
            return 0.00

    gamma_unf_SLS = Float
    def _gamma_unf_SLS_default( self ):
        if self.category == 'dead-load':
            return 1.00
        if self.category == 'additional dead-load':
            return 1.00
        if self.category == 'imposed-load':
            return 1.00


    # combination factors for imposed loads
    #
    psi_0 = Float( 1.0 )
    psi_1 = Float( 1.0 )
    psi_2 = Float( 1.0 )

    # stress resultants
    #
    def _read_state_data( self, file_name ):
        '''to read the stb-stress_resultants save the xls-worksheet 
        to a csv-file using ';' as filed delimiter and ' ' (blank)
        as text delimiter.
        '''
        print '*** read stata data ***'
        ######################################################
        # this method returns only the MAXIMUM VALUES!!!
        # @todo: dublicate the elem number and coordinates and add also the minimum values
        ######################################################

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
        # read stress_resultant csv-input file line by line in steps of two
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

        geo_columns = [ 'elem_no', 'X', 'Y', 'Z', 'D_elem' ]
        sr_columns = ['mx', 'my', 'mxy', 'nx', 'ny', 'nxy']

        return { 'elem_no' : elem_no, 'X' : X, 'Y' : Y, 'Z' : Z,
                 'mx' : mx, 'my' : my, 'mxy' : mxy,
                 'nx' : nx, 'ny' : ny, 'nxy' : nxy,
                 'geo_columns' : geo_columns, 'sr_columns' : sr_columns }

    state_data_dict = Property( Dict, depends_on = 'file_name' )
    @cached_property
    def _get_state_data_dict( self ):
        return self._read_state_data( self.file_name )

    sr_columns = List( Str )
    def _sr_columns_default( self ):
        return self.state_data_dict['sr_columns']

    sr_arr = Array
    def _sr_arr_default( self ):
        '''return the stress resultants of the loading case 
        stacked in a column array.
        '''
        sd_dict = self.state_data_dict
        return hstack( [ sd_dict[ sr_key ] for sr_key in self.sr_columns ] )

    geo_columns = List( Str )
    def _geo_columns_default( self ):
        return self.state_data_dict['geo_columns']

    geo_arr = Array
    def _geo_arr_default( self ):
        '''return the stress resultants of the loading case 
        stacked in a column array.
        '''
        sd_dict = self.state_data_dict
        return hstack( [ sd_dict[ geo_key ] for geo_key in self.geo_columns ] )


class LCManager( HasTraits ):
    '''Generate and sort the Loading case combinations
    of all specified loading cases
    '''

    #-------------------------------
    # Define loading cases:
    #-------------------------------

    lc_list = List( Instance( LC ) )
    def _lc_list_default( self ):
        return [ LC( name = 'G', category = 'dead-load', file_name = 'input_data_stress_resultants.csv' ),

                 LC( name = 'G_A', category = 'additional dead-load', file_name = 'input_data_stress_resultants.csv' ),

                 LC( name = 'Q1', category = 'imposed-load', file_name = 'input_data_stress_resultants.csv',
                     exclusive_to = ['Q2', 'Q3'], psi_0 = 0.7, psi_1 = 0.5, psi_2 = 0.3 ),

                 LC( name = 'Q2', category = 'imposed-load', file_name = 'input_data_stress_resultants.csv',
                     exclusive_to = [], psi_0 = 0.7, psi_1 = 0.5, psi_2 = 0.3 ),

                 LC( name = 'Q3', category = 'imposed-load', file_name = 'input_data_stress_resultants.csv',
                     exclusive_to = ['Q1', 'Q2'], psi_0 = 0.7, psi_1 = 0.5, psi_2 = 0.3 ),

                 LC( name = 'Q4', category = 'imposed-load', file_name = 'input_data_stress_resultants.csv',
                     exclusive_to = ['Q5'], psi_0 = 0.7, psi_1 = 0.5, psi_2 = 0.3 ),

                 LC( name = 'Q5', category = 'imposed-load', file_name = 'input_data_stress_resultants.csv',
                     exclusive_to = ['Q4'], psi_0 = 0.7, psi_1 = 0.5, psi_2 = 0.3 ),

                 LC( name = 'Q6', category = 'imposed-load', file_name = 'input_data_stress_resultants.csv',
                     exclusive_to = [], psi_0 = 0.7, psi_1 = 0.5, psi_2 = 0.3 )
                 ]

    #-------------------------------
    # check consistency
    #-------------------------------

    sr_columns = Property( List( Str ), depends_on = 'lc_list' )
    def _get_sr_columns( self ):
        '''derive the order of the stress resultants
        from the first element in 'lc_list'. The internal 
        consistency is checked separately in the 
        'check_consistency' method.
        '''
        return self.lc_list[0].sr_columns

    geo_columns = Property( List( Str ), depends_on = 'lc_list' )
    def _get_geo_columns( self ):
        '''derive the order of the stress resultants
        from the first element in 'lc_list'. The internal 
        consistency is checked separately in the 
        'check_consistency' method.
        '''
        return self.lc_list[0].geo_columns

#    @on_trait_change( 'lc_list' )
    def _check_for_consistency( self ):
        ''' check input files for consitency:
        '''
        for lc in self.lc_list[1:]:

            # check geo_columns
            #
            if lc.geo_columns != self.geo_columns:
                print 'lc.geo_columns', lc.geo_columns
                print 'self.geo_columns', self.geo_columns

                raise ValueError, 'geo_columns in loading case %s and loading case %s are not identical. Check input files for consistency!' \
                        % ( lc_list[0].name, lc.name )

            # check geo_arr
            #
            if any( lc.geo_arr ) != self.geo_arr:
                raise ValueError, 'geo_arr in loading case %s and loading case %s are not identical. Check input files for consistency!' \
                        % ( lc_list[0].name, lc.name )

            # check sr_columns
            #
            if lc.sr_columns != self.sr_columns:
                raise ValueError, 'sr_columns in loading case %s and loading case %s are not identical. Check input files for consistency!' \
                        % ( lc_list[0].name, lc.name )

            else:
                print '*** input files checked for consistency (OK) ***'
                return True

    #-------------------------------
    # lc_arr
    #-------------------------------

    lc_arr = Array
    def _lc_arr_default( self ):
        '''stack stress resultants arrays of all loading cases together.
        This yields an array of shape ( n_lc, shape, n_sr )
        '''
        sr_arr_list = [ lc.sr_arr for lc in self.lc_list ]
        return array( sr_arr_list )

    #-------------------------------
    # loading case combinations
    #-------------------------------

    def _product( self, args ):
        """
        Get all possible permutations of the security factors
        without changing the order of the loading cases.
        The method corresponds to the build-in function 'itertools.product'.
        Instead of returning a generator object a list of all 
        possible permutations is returned. As argument a list of list 
        needs to be defined. In the original version of 'itertools.product' 
        the function takes a tuple as argument ("*args"). 
        """
        pools = map( tuple, args ) #within original version args defined as *args
        result = [[]]
        for pool in pools:
            result = [x + [y] for x in result for y in pool]
        return result

    #-------------------------------
    # combi_arr - SLS
    #-------------------------------

    # Possible definitions of the serviceability limit state
    # are: 'rare', 'freq', 'perm'
    #
    combination_SLS = Str
    def _combination_SLS_default( self ):
        return 'rare'
#        return 'freq'
#        return 'perm'

    combi_arr_SLS = Property( Array )
    def _get_combi_arr_SLS( self ):
        '''array containing the security and combination factors 
        corresponding to the specified loading cases.
        This yields an array of shape ( n_lcc, n_lc )
        '''
        if self.combination_SLS != 'rare'\
            and self.combination_SLS != 'freq'\
            and self.combination_SLS != 'perm':
            raise ValueError, "SLS - combination % s does not exist ! " % ( self.combination_SLS )
        else:
            print " *** combination of SLS for % s *** " % ( self.combination_SLS )

        #-------------------------------------------------------------
        # 0. step: get auxillary lists for permutations:
        #          ('category_list', 'psi_list', 'imposed_idx_list', 'gamma_list')
        #-------------------------------------------------------------

        # list containing the type of load for each loading case defined in 'lc_list'
        #
        category_list = []

        # list containing the value of 'psi_0' for each loading case defined in 'lc_list'
        # (for body-loads a value of 1.0 is added)
        #
        psi_list = []

        # list containing indices of the imposed loading cases
        # (order taken from 'lc_list')
        #
        imposed_idx_list = []

        # definition of load cases and values of psi for SLS
        #
        if self.combination_SLS == 'rare':
            # 'rare' combination defined as:
            # G + Q + sum( psi_0i * Qi )         
            #
            for i_lc, lc in enumerate( self.lc_list ):
                cat = lc.category
                if cat == 'imposed-load':
                    category_list.append( 'imposed-load' )
                    psi_list.append( lc.psi_0 )
                    imposed_idx_list.append( i_lc )
                elif cat == 'dead-load' or lc.category == 'additional dead-load':
                    category_list.append( 'dead-load' )
                    psi_list.append( 1.0 )
                else:
                    raise ValueError, "'category' of loading case instance % s is undefined or does not exist ! " % ( lc.name )

        if self.combination_SLS == 'freq' or self.combination_SLS == 'perm':
            # 'frequent' combination defined as:
            # G + psi_1 * Q + sum( psi_2i * Qi )         
            #
            for i_lc, lc in enumerate( self.lc_list ):
                cat = lc.category
                if cat == 'imposed-load':
                    category_list.append( 'imposed-load' )
                    psi_list.append( lc.psi_2 )
                    imposed_idx_list.append( i_lc )
                elif cat == 'dead-load' or lc.category == 'additional dead-load':
                    category_list.append( 'dead-load' )
                    psi_list.append( 1.0 )
                else:
                    raise ValueError, "'category' of loading case instance % s is undefined or does not exist ! " % ( lc.name )

        # for the SLS as combination the gamma factors do not exist.
        # They are used here, and only hold for 1 or 0, to distinguish
        # favorable or unfavorable imposed load combinations.
        #
        gamma_list = []
        for i_lc, lc in enumerate( self.lc_list ):
            if category_list[ i_lc ] == 'imposed-load':
                gamma_list.append( [ lc.gamma_fav_SLS, lc.gamma_unf_SLS ] )
            else:
                gamma_list.append( [ 1.0 ] )

        #--------------------------------------------------------
        # 1. step: get all permutations of the security factors:
        #          ('gamma_factors_SLS')
        #--------------------------------------------------------

        # permutations of all gamma factors
        #
        permutation_list = self._product( gamma_list )

        # combination array without psi values and exclusive load cases
        #
        combi_arr = array( permutation_list )

        #--------------------------------------------------------
        # 2. step: permutate the leading imposed load and multiply the
        #          remaining imposed-loads with combination factors ('psi')
        #--------------------------------------------------------

        # if no imposed-loads have been defined return the combination array
        # as obtained in step 1:
        # 
        if imposed_idx_list == []:
            return combi_arr

        # for each imposed load case, multiply all other with psi value (non leading) 
        # and evaluate comb_arr for leading imposed load, where the leading imposed load
        # is unequal to zero
        # 
        combi_arr_psi_list = []
        if self.combination_SLS == "perm":
            psi_arr = array( psi_list )
            combi_arr_psi = combi_arr * transpose( psi_arr )

        elif self.combination_SLS != "perm":
            for idx in imposed_idx_list:
                psi_arr_i = array( psi_list )
                # set imposed load case i as lead
                if self.combination_SLS == "rare":
                    # define leading value  = 1.0
                    #
                    psi_arr_i[ idx ] = 1.0
                if self.combination_SLS == "freq":
                    # define leading value = psi_1
                    #
                    psi_arr_i[ idx ] = self.lc_list[ idx ].psi_1

                psi_arr_i = transpose( psi_arr_i )
                combi_arr_lead_i = combi_arr[where( combi_arr[ :, idx ] != 0 )] * psi_arr_i
                combi_arr_psi_list.append( combi_arr_lead_i )

            combi_arr_psi_no_0 = vstack( combi_arr_psi_list )

            # missing cases without any dead load have to be added 
            # get combinations with all imposed = 0 
            #
            lcc_all_imposed_zero = where( ( combi_arr[:, imposed_idx_list] == 0 ).all( axis = 1 ) )

            # add to combinations
            #
            combi_arr_psi = vstack( ( combi_arr[lcc_all_imposed_zero], combi_arr_psi_no_0 ) )

        #--------------------------------------------------------
        # 3. step: delete the combinations containing more then
        #          one exclusive load case ('exclusive_to):
        #--------------------------------------------------------

        # get exclusive values , which are related to each other
        #
        exclusive_list = []
        for i in range( 0, self.n_lc ):
            # check for exclusive defined
            #
            if self.lc_list[i].exclusive_to != []:
                # get related load case number
                #
                for j in range ( 0, self.n_lc ):
                    # multiple exclusive load cases for one load case have to be considered
                    #
                    for exclusive_i in self.lc_list[i].exclusive_to:
                        if exclusive_i == self.lc_list[j].name:
                            exclusive_list.append( sorted( [i, j] ) )

        # remove entries if they exist more then once
        # 
        exclusive_list_unique = []
        for x in exclusive_list:
            if x not in exclusive_list_unique:
                 exclusive_list_unique.append( x )

        # remove combinations that are defined exclusive
        # 
        combi_arr_psi_exclusive = combi_arr_psi
        for exclusive_case in exclusive_list_unique:
            # check where maximum one value of the exclusive load cases is unequal to one 
            # f.e.         1.5 0.9 0.8 all exclusive
            #        ->    1.0 1.0 1.0 -> sum = 3
            #        ->    false -> non accepted combination
            # f.e.         1.5 0.0 0.0 all exclusive
            #        ->    1.0 0.0 0.0 -> sum = 1
            #        ->    true -> accepted combination
            #
            true_combi = where( sum( combi_arr_psi_exclusive
                                   [:, exclusive_case] != 0, axis = 1.0 ) <= 1.0 )
            combi_arr_psi_exclusive = combi_arr_psi_exclusive[ true_combi ]

        return combi_arr_psi_exclusive

    #-------------------------------
    # combi_arr - ULS
    #-------------------------------

    combi_arr_ULS = Property( Array )
    def _get_combi_arr_ULS( self ):
        '''return an array containing the security and combination
        factors corresponding to the specified loading cases.
        This yields an array of shape ( n_lcc, n_lc ).
        '''
        #-------------------------------------------------------------
        # 0. step: get auxillary lists for permutations:
        #          ('category_list', 'psi_list', 'imposed_idx_list', 'gamma_list')
        #-------------------------------------------------------------

        # list containing the type of load for each loading case defined in 'lc_list'
        #
        category_list = []

        # list containing the value of 'psi_0' for each loading case defined in 'lc_list'
        # (for body-loads a value of 1.0 is added)
        #
        psi_list = []

        # list containing indices of the imposed loading cases
        # (order taken from 'lc_list')
        #
        imposed_idx_list = []

        for i_lc, lc in enumerate( self.lc_list ):
            cat = lc.category
            if cat == 'imposed-load':
                category_list.append( 'imposed-load' )
                psi_list.append( lc.psi_0 )
                imposed_idx_list.append( i_lc )
            elif cat == 'dead-load' or lc.category == 'additional dead-load':
                category_list.append( 'dead-load' )
                psi_list.append( 1.0 )
            else:
                raise ValueError, "'category' of loading case instance % s is undefined or does not exist ! " % ( lc.name )

        # list of list containing gamma factors for all loading cases:
        #
        gamma_list = [ [ lc.gamma_fav, lc.gamma_unf ] for lc in self.lc_list ]

        #--------------------------------------------------------
        # 1. step: get all permutations of the security factors:
        #          ('gamma_factors')
        #--------------------------------------------------------

        # permutations of all gamma factors
        #
        permutation_list = self._product( gamma_list )

        # combination array without psi values and exclusive load cases
        #
        combi_arr = array( permutation_list )

        #--------------------------------------------------------
        # 2. step: permutate the leading imposed load and multiply the
        #          remaining imposed-loads with combination factors ('psi0')
        #--------------------------------------------------------

        # if no imposed-loads have been defined return the combination array
        # as obtained in step 1:
        # 
        if imposed_idx_list == []:
            return combi_arr

        # for each imposed load case, multiply all other lc's with their psi value (non leading) 
        # and evaluate 'combi_arr' for leading imposed load, where the leading imposed load
        # is unequal to zero
        # 
        combi_arr_psi_list = []
        for idx in imposed_idx_list:
            psi_arr_i = array( psi_list )
            # set imposed load case i_lc as lead
            psi_arr_i[ idx ] = 1.0
            psi_arr_i = transpose( psi_arr_i )
            combi_arr_lead_i = combi_arr[ where( combi_arr[ :, idx ] != 0 ) ] * psi_arr_i
            combi_arr_psi_list.append( combi_arr_lead_i )

        combi_arr_psi_no_0 = vstack( combi_arr_psi_list )

        # missing cases without any imposed load have to be added 
        # get combinations with all imposed set to 0 
        #
        lcc_all_imposed_zero = where( ( combi_arr[:, imposed_idx_list] == 0 ).all( axis = 1 ) )

        # add to combinations
        #
        combi_arr_psi = vstack( ( combi_arr[ lcc_all_imposed_zero ], combi_arr_psi_no_0 ) )

        #--------------------------------------------------------
        # 3. step: delete the combinations containing more then
        #          one exclusive load case ('exclusive_to):
        #--------------------------------------------------------

        # get exclusive values , which are related to each other
        #
        exclusive_list = []
        for i in range( 0, self.n_lc ):
            # check for exclusive defined
            #
            if self.lc_list[i].exclusive_to != []:
                # get related load case number
                #
                for j in range ( 0, self.n_lc ):
                    # multiple exclusive load cases for one load case have to be considered
                    #
                    for exclusive_i in self.lc_list[i].exclusive_to:
                        if exclusive_i == self.lc_list[j].name:
                            exclusive_list.append( sorted( [i, j] ) )

        # remove entries if they exist more then once
        # 
        exclusive_list_unique = []
        for x in exclusive_list:
            if x not in exclusive_list_unique:
                 exclusive_list_unique.append( x )

        # remove combinations that are defined exclusive
        # 
        combi_arr_psi_exclusive = combi_arr_psi
        for exclusive_case in exclusive_list_unique:
            # check where maximum one value of the exclusive load cases is unequal to one 
            # f.e.         1.5 0.9 0.8 all exclusive
            #        ->    1.0 1.0 1.0 -> sum = 3
            #        ->    false -> non accepted combination
            # f.e.         1.5 0.0 0.0 all exclusive
            #        ->    1.0 0.0 0.0 -> sum = 1
            #        ->    true -> accepted combination
            #
            true_combi = where( sum( combi_arr_psi_exclusive
                                   [:, exclusive_case] != 0, axis = 1.0 ) <= 1.0 )
            combi_arr_psi_exclusive = combi_arr_psi_exclusive[true_combi]

        return combi_arr_psi_exclusive

    #-------------------------------
    # Array dimensions:
    #-------------------------------

    n_sr = Property( Int )
    def _get_n_sr( self ):
        return len( self.sr_columns )

    n_lc = Property( Int )
    def _get_n_lc( self ):
        return len( self.lc_list )

    n_lcc = Property( Int )
    def _get_n_lcc( self ):
        return self.combi_arr_ULS.shape[0]

    shape = Property( Int )
    def _get_shape( self ):
        return self.lc_list[0].sr_arr.shape[0]

    lcc_arr_ULS = Property( Array )
    def _get_lcc_arr_ULS( self ):
        """
        Array of all ULS loading case combinations
        This yields an array of shape ( n_lcc, shape, n_sr )
        """
        return self._get_lcc_arr( self.combi_arr_ULS )

    #-------------------------------
    # lcc_arr - SLS
    #-------------------------------

    lcc_arr_SLS = Property( Array )
    def _get_lcc_arr_SLS( self ):
        """
        Array of all SLS loading case combinations
        This yields an array of shape ( n_lcc, shape, n_sr )
        """
        return self._get_lcc_arr( self.combi_arr_SLS )

    #-------------------------------
    # lcc_arr - ULS
    #-------------------------------

    def _get_lcc_arr( self, combi_arr ):
        '''Array of all loading case combinations following the
        loading cases define in 'lc_list' and the combinations
        defined in 'combi_arr'.
        This yields an array of shape ( n_lcc, shape, n_sr ) 
        '''
        # 'combi_arr' is of shape ( n_lcc, n_lc )
        # 'lc_arr' is of shape ( n_lc, shape, n_sr )
        #
        lc_arr = self.lc_arr

        # Broadcasting is used to generate array containing the multiplied lc's
        # yielding an array of shape ( n_lcc, n_lc, shape, n_sr )
        #
        lc_combi_arr = lc_arr[ newaxis, :, :, :] * combi_arr[:, :, newaxis, newaxis]

        # Then the sum over index 'n_lc' is evaluated yielding 
        # an array of all loading case combinations.
        # This yields an array of shape ( n_lcc, shape, n_sr ) 
        #
        lcc_arr = sum( lc_combi_arr, axis = 1.0 )

        # create an array for identification of the load case combination
        # This yields an array of shape ( n_lcc, shape, 1 ) 
        #
        n_lcc = self.n_lcc
        shape = self.shape
        idx_col = ones( ( n_lcc, shape, 1 ) ) * arange( n_lcc )[ :, newaxis, newaxis ]

        # Note that 'vstack' must have identical shape in the first dimension, therefore 
        # 'transpose' is used to reverse the order of 'lcc' and 'idx_col'. After 'vstack'
        # has been performed the order is reversed back to the former order.
        #
        lcc = transpose( vstack( ( transpose( lcc_arr ), transpose( idx_col ) ) ) )

        return lcc


    #-------------------------------
    # min/max-values - for verification only!
    #-------------------------------

    def get_min_max( self, lcc_arr_i ):
        """ get the surrounding curve of all 'lcc' values
        """
#        lcc_arr = self.lcc_arr

        arg_min = argmin( lcc_arr_i[:, :, :-1], axis = 0 )
        min_val = min( lcc_arr_i[:, :, :-1], axis = 0 )
        arg_max = argmax( lcc_arr_i[:, :, :-1], axis = 0 )
        max_val = max( lcc_arr_i[:, :, :-1], axis = 0 )

        n_elemens = shape( lcc_arr_i )[1]

        min_max_arr = vstack( ( min_val, arg_min, max_val, arg_max ) )
        min_max_arr = min_max_arr.reshape( shape( min_max_arr )[0]
                                          / n_elemens, n_elemens, 6 )

        return min_max_arr

    #-------------------------------
    # lcc_tree 
    #-------------------------------

    lcc_tree = Dict
    def _lcc_tree_default( self ):
        '''dictionary containing the loading case combinations
        for the different limit states.
        '''
        sr_columns = self.sr_columns

        # return a dictionary of the stress resultants
        # this is used by InfoShellPro to determine the stress 
        # resultants of the current limit state 
        #
        lcc_columns = sr_columns + ['combi_key']

        lcc_dict = {}
        for ls in ['ULS', 'SLS']:

            if ls == 'ULS':
                lcc_arr = self.lcc_arr_ULS
            if ls == 'SLS':
                lcc_arr = self.lcc_arr_SLS

            # this yields a 2D - array of shape (n_lcc * shape, n_sr )
            lcc_arr_2D = lcc_arr.reshape( self.n_lcc * self.shape, self.n_sr + 1 )

            # this is equivalent to: 
            # lcc_arr_2D = vstack( [ lcc_arr[i_lcc, :, :] for i_lcc in range( self.n_lcc )] )

            col_dict = {}
            for idx, name in enumerate( lcc_columns ):
                col_dict[ name ] = lcc_arr_2D[ :, idx ]
            lcc_dict[ ls ] = col_dict

        return lcc_dict

if __name__ == '__main__':
    lcm = LCManager()
    lcm.combination_SLS = 'freq'
    print shape( lcm.combi_arr_SLS )

#    print 'combi_arr.shape', lcm.combi_arr_ULS.shape
    print 'combi_arr_SLS', lcm.combi_arr_SLS
    print 'combi_arr_ULS', lcm.combi_arr_ULS
#    print 'combi_arr', lcm.combi_arr_ULS

#    print 'combi_arr_SLS', lcm.combi_arr_ULS

#    print 'lc_arr.shape', lcm.lc_arr.shape
#    print lcm.lc_list[0].sr_arr
#    print 'sr_arr.shape', lcm.lc_list[0].sr_arr.shape
#    print 'lcc_arr_ULS.shape', lcm.lcc_arr_ULS.shape
#    print 'lcc_tree[ULS].shape', lcm.lcc_tree['SLS']['mx']





