


from enthought.traits.api import HasTraits, Float, Property, cached_property, \
    Event, Array, Instance, File, Int, Directory, Button, Range, Enum, \
     on_trait_change, Bool, Trait, HasPrivateTraits, Constant, List, Tuple, \
     DelegatesTo
from enthought.traits.ui.api import Item, View, Group, Handler, HSplit, VSplit, \
    HGroup, HSplit, VGroup, Tabbed, Label, Controller, ModelView, \
    FileEditor, DirectoryEditor, \
    HistoryEditor, RangeEditor
from enthought.traits.ui.menu import OKButton, CancelButton, Action, Menu, \
    MenuBar
from enthought.tvtk.pyface.scene_editor import SceneEditor
from matplotlib.figure import Figure
from numpy import loadtxt, min, array, mean, std, arange, histogram, zeros_like, \
    meshgrid, ones_like, histogram2d, c_, r_, equal, cumsum, vstack, hstack, savetxt, \
    sqrt, sum, all, zeros_like, zeros, ones, diff, where, unique, isnan, pi, invert
from numpy.random import random
from os.path import join
from promod.simdb import SimDB
from traits.editors.mpl_figure_editor import MPLFigureEditor
import matplotlib.pyplot as plt
import numpy.ma as ma
import os
import re
from enthought.traits.trait_types import DelegatesTo
from time import time
from scipy.sparse import csr_matrix

#import wxversion

#wxversion.select( '2.8' )

simdb = SimDB()
data_dir = join( simdb.exdata_dir, 'trc', 'yarn_structure', 'VET', 'raw_data' )
#data_dir = join( simdb.exdata_dir, 'trc', 'yarn_structure', 'MAG', 'raw_data' )
#data_dir = join( simdb.exdata_dir, 'trc', 'yarn_structure', 'TEST' )


class YMBRawData( HasTraits ):
    '''
        Read raw data files and convert length units px to micrometers
    '''
    raw_data_directory = Directory( data_dir, entries=10,
                                    changed_source=True,
                                    auto_set=False, enter_set=True )

    # constant -- convert pixels to micrometers
    px2mum = Constant( 1.2, input=True )

    cut_raw_data = Property( Tuple, depends_on='raw_data_directory' )
    @cached_property
    def _get_cut_raw_data( self ):
        '''
            Given a directory name, returns Tuple of NumPy arrays of variables 
            that can be obtained from raw data (y, z, r, d, contact fraction).
            Cut filenames must satisfy the following format '*Schnitt(num).txt'
            Return:
                y,z -- 1D arrays of coordinates
                r -- 1D array of diameters 
                d -- 1D array of distances from the edge
                cf -- 2D array of contact fraction components
                offset -- 1D array of number of filaments in the cut
        '''
        files = []
        num = []
        paths = os.listdir( self.raw_data_directory )  # list of paths in that dir
        for fname in paths:
            match = re.search( r'\w+Schnitt(\d+).txt', fname ) # find files of given pattern
            if match:
                files.append( os.path.abspath( os.path.join( self.raw_data_directory, fname ) ) )
                num.append( int( match.group( 1 ) ) )
        num, files = zip( *sorted( zip( num, files ) ) )
        y_arr = array( [] )
        z_arr = array( [] )
        r_arr = array( [] )
        d_arr = array( [] )
        cf_component_arr = array( [] ).reshape( 0, 36 )
        offset = array( [] )

        for n, file in zip( num, files ):
            print 'reading cut data_file -- slice %i' % n
            data = loadtxt( file ,
                              skiprows=1,
                              usecols=None )

            y = data[ :, 1] / self.px2mum
            z = data[ :, 2] / self.px2mum
            r = data[ :, 3] / self.px2mum
            d = data[ :, 4] / self.px2mum
            cf_component = data[ :, 5:41]

            y_arr = hstack( [y_arr, y] )
            z_arr = hstack( [ z_arr, z] )
            r_arr = hstack( [ r_arr, r] )
            d_arr = hstack( [ d_arr, d] )
            cf_component_arr = vstack( [cf_component_arr, cf_component] )
            offset = hstack( [offset, len( y )] )

        y_shift = min( y_arr )
        z_shift = min( z_arr )
        if y_shift < 0:
            print 'y-coordinites shifted by %f' % y_shift
            y_arr = y_arr - y_shift
        if z_shift < 0:
            print 'z-coordinites shifted by %f' % z_shift
            z_arr = z_arr - z_shift
        return y_arr, z_arr, r_arr, d_arr, cf_component_arr, cumsum( offset )

    # number of cuts in the raw data (numbering from zero)
    n_cuts = Property( Int, depends_on='+changed_source' )
    @cached_property
    def _get_n_cuts( self ):
        n_cuts = self.connectivity.shape[1] - 1
        print 'number of cuts', n_cuts
        return n_cuts

    # number of cuts in the raw data (numbering from zero)
    n_filaments = Property( Int, depends_on='+changed_source' )
    @cached_property
    def _get_n_filaments( self ):
        #y_list = self.cut_raw_data[0]
        n_filaments = self.connectivity.shape[0] - 1 # max( [ cut.shape[0] for cut in y_list ] )
        print 'number of filaments', n_filaments
        return n_filaments

    # yarn cut coordinates
    # TODO: change with new data sets? probably yes
    cut_x = Property( depends_on='+changed_source' )
    @cached_property
    def _get_cut_x( self ):
        '''
            Read data with cut x-coordinates
            First cut has x-coordinate = 0
            Return NumPy array
        '''
        print 'reading data_file -- slice_height.txt'
        #slice_height_file = join( self.raw_data_directory, 'slice_height.txt' )
        cut_x_file = join( self.raw_data_directory, 'cut_x.txt' )

        cut_x = loadtxt( cut_x_file, dtype='float',
                                skiprows=1, usecols=( 1, 3 ) ) / self.px2mum
        cut_x = unique( cut_x.flatten() )
        print 'cut coordinates', cut_x
        return cut_x
        #return cumsum( slice_height[:-1] )
        #return cumsum( array( [500] * 15 ) ) / self.px2mum

    # filament connections between cuts
    connectivity = Property( Array, depends_on='+changed_source' )
    @cached_property
    def _get_connectivity( self ):
        '''
            Read connectivity data.
            Return NumPy array
        '''
        print 'reading data_file -- connectivity.txt'
        connectivity_file = join( self.raw_data_directory, 'connectivity.txt' )
        connect = loadtxt( connectivity_file, dtype='float', skiprows=1 )
        # replace 0 by -1 (not identified)
        # TODO: will be changed in raw data
        connect[ equal( 0, connect )] = -1
        return connect

    #traits_view = View( Item( 'raw_data_directory' ) )


class YMBData( YMBRawData ):
    '''
        Preprocessing raw data
    '''
    #raw_data = Instance( YarnRawData, () )

    #cut_raw_data = DelegatesTo( 'raw_data' )
    #slice_height = DelegatesTo( 'raw_data' )
    #connectivity = DelegatesTo( 'raw_data' )
    input_change = Event
    @on_trait_change( '+changed_source, +changed_config' )
    def _set_input_change( self ):
        print '*** raising input change in CTT'
        self.input_change = True

    cf_limit = Range( value=.2, low=0.0,
                      changed_config=True,
                      high=1.0, enter_set=True, auto_set=False )

    cut_data = Property( Tuple, depends_on='+changed_source' )
    @cached_property
    def _get_cut_data( self ):
        '''
            Prepare grid from raw data
            Output: 2D NumPy arrays (x, y, z, r, d, cf) 
            (cut (columns) x filament (rows))
        '''
        print 'preparing cut_data'
        y_arr, z_arr, r_arr, d_arr, cf_arr, offset = self.cut_raw_data
        slice_height_arr = self.cut_x
        filam_connect_arr = self.connectivity

        offset_arr = zeros_like( offset )
        offset_arr[1:] = offset[:-1]

        offset_arr = ones_like( filam_connect_arr ) * offset_arr
        map_arr = zeros_like( filam_connect_arr )
        map_arr[invert( self.mask_arr )] = filam_connect_arr[invert( self.mask_arr )] + offset_arr[invert( self.mask_arr )]
        map_arr[self.mask_arr] = filam_connect_arr[self.mask_arr] * 0 - 1
        map_arr = array( map_arr, dtype='int' )

        # stack 1d arrays in list into single array with -1 at the end
        #x_raw_arr = hstack( [hstack( ones_like( map_arr ) * slice_height_arr ), [-1]] )
        y_arr = hstack( [y_arr, [-1]] )
        z_arr = hstack( [z_arr, [-1]] )
        r_arr = hstack( [r_arr, [-1]] )
        d_arr = hstack( [d_arr, [-1]] )

        # solve cf from contact fraction components (0-contact with matrix,
        # 1-contact with fiber, 2-no contact) -- 36 values around fiber diameter
        cf_arr = sum( cf_arr == 0, axis=1 ) / 36.
        cf_arr = hstack( [ cf_arr , [-1]] )

        # map data arrays according to connection data
        x_arr = ones_like( map_arr ) * slice_height_arr
        x_arr[self.mask_arr] = x_arr[self.mask_arr] * 0 - 1
        y_arr = y_arr[map_arr]
        z_arr = z_arr[map_arr]
        r_arr = r_arr[map_arr]
        d_arr = d_arr[map_arr]
        cf_arr = cf_arr[map_arr]

        print 'array shapes', x_arr.shape, y_arr.shape, z_arr.shape, r_arr.shape, d_arr.shape, cf_arr.shape
        return x_arr, y_arr, z_arr, r_arr, d_arr, cf_arr

    mask_arr = Property( Array, depends_on='+changed_source' )
    @cached_property
    def _get_mask_arr( self ):
        print 'preparing mask_array'
        return self.connectivity < 0

    radius = Property( Array, depends_on='+changed_source' )
    @cached_property
    def _get_radius( self ):
        print 'preparing radius'
        return self.cut_data[3]

    edge_distance = Property( Array, depends_on='+changed_source' )
    @cached_property
    def _get_edge_distance( self ):
        print 'preparing edge_distance'
        return self.cut_data[4]

    cf = Property( Array, depends_on='+changed_source' )
    @cached_property
    def _get_cf( self ):
        print 'preparing contact fraction'
        #savetxt( 'cf.txt', self.cut_data[5] )
        return self.cut_data[5]

    # filament area
    area = Property( Array, depends_on='+changed_source' )
    @cached_property
    def _get_area( self ):
        area = pi * self.radius ** 2
        area[self.length_between_cuts.mask] = -1
        return area

    # yarn area in the cut
    cut_area = Property( Array, depends_on='+changed_source' )
    @cached_property
    def _get_cut_area( self ):
        area = pi * self.radius ** 2
        area[self.length_between_cuts.mask] = 0
        area = sum( area, axis=0 )
        return area.reshape( 1, self.n_cuts + 1 )

    # @kelidas: Old slack for whole fiber length
#    slack = Property( depends_on='cut_data, length_between_cuts' )
#    @cached_property
#    def _get_slack( self ):
#        '''
#            return slack of the whole filament
#        '''
#        mask = self.length_between_cuts.mask
#        length_between_cuts = array( self.length_between_cuts )
#        length_between_cuts[mask] = 0
#        fib_len = sum( length_between_cuts, axis=1 )
#        cut_len = zeros_like( self.slice_height )
#        cut_len[1:] = cut_len[1:] + diff( self.slice_height )
#        c_dist = ones_like( length_between_cuts ) * cut_len[1:]
#        c_dist[mask] = 0
#        c_len = sum( c_dist, axis=1 )
#        slack = ( ( fib_len - c_len ) / fib_len )
#        slack = slack.reshape( len( slack ), 1 )
#        return ones_like( self.cf ) * array( slack )

    # real free filament length (sum of the lengths without matrix contact
    # obtained by linear interpolation between cuts)
    free_length_real = Property( Array, depends_on='+changed_source, +changed_config' )
    @cached_property
    def _get_free_length_real( self ):
        start = time()
        print 'preparing free_length_real'
        cf = self.cf.copy()
        cf[cf <= self.cf_limit] = 0
        #cfl = ma.array( cfl, mask=self.length_between_cuts.mask )
        cf[cf != 0] = 1 # contact
        cf[cf == 0] = 0 #  no contact

        cfl = cf[:, :-1] + cf[:, 1:]
        cfl[cfl == 0] = -1
        cfl[cfl > 0] = 1
        #savetxt( 'temp.txt', cfl )
        cfl_id = cfl[:, :-1] * cfl[:, 1:]
        #print cfl_id
        # length between cuts
        #cfl = cfl * diff( self.cut_x )
        cfl = cfl * array( self.length_between_cuts )
        cfl[self.length_between_cuts.mask] = 0
        #print cfl
        # TODO: possible--use numpy along_axis
        cut_pos = range( 0, self.n_cuts )
        cfl_cuts = []
        for cut in cut_pos:
            cfl_cut = []
            for i, il in zip( cfl_id, cfl ):
                switches = where( i < 0 )[0]
                switches = hstack( [array( [-1] ), switches, [len( il ) - 1]] )
                #try:
                idx = where( switches >= cut )[0][0]
                right_limit = switches[ idx ] + 1
                left_limit = switches[ idx - 1 ] + 1
                #print left_limit, right_limit
                #print left_limit, right_limit
                #print il[ left_limit: right_limit ]
                if len( il[ left_limit: right_limit ] ) != 0 and all( il[ left_limit: right_limit ] <= 0 ) == True:
                    cfl_cut.append( sum( abs( il[ left_limit: right_limit ] ) ) )
                else:
                    cfl_cut.append( 0 )
                #except:
                #    cfl_cut.append( 0 )
            #print cfl_cut
            #print len( array( cfl_cut ) )
            cfl_cuts.append( array( cfl_cut ) )
        #print vstack( cfl_cuts ).T
        cfl_cuts = vstack( cfl_cuts ).T
        cfl_cuts[self.length_between_cuts.mask] = -1
        #print cfl_cuts
        print time() - start
        print cfl_cuts
        print cfl_cuts.shape
        return cfl_cuts

# cut free filament length (sum of the cut distances without matrix contact)
    free_length_cut = Property( Array, depends_on='+changed_source, +changed_config' )
    @cached_property
    def _get_free_length_cut( self ):
        start = time()
        print 'preparing free_length_cut'
        cf = self.cf.copy()
        cf[cf <= self.cf_limit] = 0
        #cfl = ma.array( cfl, mask=self.length_between_cuts.mask )
        cf[cf != 0] = 1 # contact
        cf[cf == 0] = 0 #  no contact

        cfl = cf[:, :-1] + cf[:, 1:]
        cfl[cfl > 0] = -1
        cfl[cfl == 0] = 1
        
        data_bas = cfl
        x_data = self.cut_x # linspace( 0, 10, data_bas.shape[1] )

        data = hstack( [ -ones( ( data_bas.shape[ 0 ], 1 ) ), data_bas / abs( data_bas ),
                         - ones( ( data_bas.shape[ 0 ], 1 ) ) ] )
        print 'data', data
    
        r = range( 0, data.shape[0] )
        filaments, switches = where( data[:, 0:-1] * data[:, 1:] < 0 )
        print 'filaments', filaments
        print 'switches', switches
        n_segments = len( filaments ) / 2.
        fil_idx_arr = filaments.reshape( ( n_segments, 2 ) )
        switch_idx_arr = switches.reshape( ( n_segments, 2 ) )
    
        print fil_idx_arr
        print switch_idx_arr
    
        segment_lengths = x_data[ switch_idx_arr[:, 1] ] - x_data[ switch_idx_arr[:, 0] ]
        print 'segment lengths', segment_lengths
    
        # data for sparse matrix CSR
        data_row = segment_lengths.repeat( switch_idx_arr[:, 1] - switch_idx_arr[:, 0] )
        # row info for data
        row = fil_idx_arr[:, 0].repeat( switch_idx_arr[:, 1] - switch_idx_arr[:, 0] )
    
        print 'data_row', data_row
        print 'row', row
    
        # position assuming that data.flatten()
        switch_idx_arr = ( switch_idx_arr + fil_idx_arr * data.shape[1] ).flatten()
        print 'switch idx array', switch_idx_arr
        switch_idx_arr2 = switch_idx_arr[1:] - switch_idx_arr[:-1]
        print 'asdfdas', switch_idx_arr2
    
    
        # create array  [1,0,1,0,......]
        aran = arange( 0, len( switch_idx_arr2 ) )
        mask = aran % 2 == 0
        aran[mask] = True
        aran[invert( mask )] = False
        # repeat values
        aran = ( aran.repeat( switch_idx_arr2 ) ).astype( bool )
        print aran
    
        a = arange( min( switch_idx_arr ), max( switch_idx_arr ) )
        print a
        print 'fas', a[aran]
        col = a[aran] - row * data.shape[1]
    
        cfl_cuts = array( csr_matrix( ( data_row, ( row, col ) ), shape=data_bas.shape ).todense() )
        print cfl_cuts
        print data_bas
        
        print time() - start
        cfl_cuts[self.length_between_cuts.mask] = -1
        return cfl_cuts


#    # cut free filament length (sum of the lengths without matrix contact
#    # obtained by linear interpolation between cuts)
#    free_length_real = Property( Array, depends_on='+changed_source, +changed_config' )
#    @cached_property
#    def _get_free_length_real( self ):
#        start = time()
#        print 'preparing free_length_real'
#        cf = self.cf.copy()
#        cf[cf <= self.cf_limit] = 0
#        #cfl = ma.array( cfl, mask=self.length_between_cuts.mask )
#        cf[cf != 0] = 1 # contact
#        cf[cf == 0] = 0 #  no contact
#
#        cfl = cf[:, :-1] + cf[:, 1:]
#        cfl[cfl == 0] = -1
#        cfl[cfl > 0] = 1
#        #savetxt( 'temp.txt', cfl )
#        cfl_id = cfl[:, :-1] * cfl[:, 1:]
#        #print cfl_id
#        cfl = cfl * diff( self.cut_x )#* array( self.length_between_cuts )
#        cfl[self.length_between_cuts.mask] = 0
#        #print cfl
#        # TODO: possible--use numpy along_axis
#        cut_pos = range( 0, self.n_cuts )
#        cfl_cuts = []
#        for cut in cut_pos:
#            cfl_cut = []
#            for i, il in zip( cfl_id, cfl ):
#                switches = where( i < 0 )[0]
#                switches = hstack( [array( [-1] ), switches, [len( il ) - 1]] )
#                #try:
#                idx = where( switches >= cut )[0][0]
#                right_limit = switches[ idx ] + 1
#                left_limit = switches[ idx - 1 ] + 1
#                #print left_limit, right_limit
#                #print left_limit, right_limit
#                #print il[ left_limit: right_limit ]
#                if len( il[ left_limit: right_limit ] ) != 0 and all( il[ left_limit: right_limit ] <= 0 ) == True:
#                    cfl_cut.append( sum( abs( il[ left_limit: right_limit ] ) ) )
#                else:
#                    cfl_cut.append( 0 )
#                #except:
#                #    cfl_cut.append( 0 )
#            #print cfl_cut
#            #print len( array( cfl_cut ) )
#            cfl_cuts.append( array( cfl_cut ) )
#        #print vstack( cfl_cuts ).T
#        cfl_cuts = vstack( cfl_cuts ).T
#        cfl_cuts[self.length_between_cuts.mask] = -1
#        #print cfl_cuts
#        print time() - start
#        print cfl_cuts
#        return cfl_cuts

    slack = Property( Array, depends_on='+changed_source, +changed_config' )
    @cached_property
    def _get_slack( self ):
        '''
            return slack in the cut of the free length part of the filament
        '''
        start = time()
        print 'preparing slack'
        slack = ( self.free_length_real - self.free_length_cut ) / self.free_length_real
        slack[isnan( slack )] = 0
        slack[self.length_between_cuts.mask] = -1
        #savetxt( 'slack.txt', slack )
        print time() - start
        return slack

    # length obtained by linear interpolation between cuts
    length_between_cuts = Property( Array, depends_on='+changed_source, +changed_config' )
    @cached_property
    def _get_length_between_cuts( self ):
        '''
            Return linear filament length between cuts.
            Matrix contain mask information.
        '''
        print 'preparing length_between_cuts'
        mask = self.mask_arr
        x_arr = ma.array( self.cut_data[0], mask=mask )
        y_arr = ma.array( self.cut_data[1], mask=mask )
        z_arr = ma.array( self.cut_data[2], mask=mask )
        return sqrt( ( x_arr[:, 1:] - x_arr[:, :-1] ) ** 2
                + ( y_arr[:, 1:] - y_arr[:, :-1] ) ** 2
                + ( z_arr[:, 1:] - z_arr[:, :-1] ) ** 2 )

    traits_view = View( HGroup( Item( 'raw_data_directory', springy=True ),
                                Item( 'cf_limit' ) ) )


class YMBSlider( HasTraits ):
    '''
        Slicing data arrays (2d NumPy array)
        Return data for statistics (1d NumPy array)
    '''
    data = Instance( YMBData, () )

    n_cuts = DelegatesTo( 'data' )
    n_filaments = DelegatesTo( 'data' )

    var_enum = Trait( 'radius',
                      {'radius' : 'radius',
                       'area':'area',
                       'yarn area in the cut':'cut_area',
                       'shortest distance from the edge' : 'edge_distance',
                       'contact fraction' : 'cf',
                       'slack' : 'slack',
                       'real free length': 'free_length_real',
                       'cut free length': 'free_length_cut'} )

    zero = Constant( 0 )
    cut_slider = Int( changed_range=True )
    cut_slider_on = Bool( True, changed_range=True )

    filament_slider = Int( changed_range=True )
    filament_slider_on = Bool( False, changed_range=True )

    stat_data = Property( Array,
                          depends_on='var_enum, +changed_range, data.+changed_source, data.+changed_config' )
    @cached_property
    def _get_stat_data( self ):
        if self.cut_slider_on == True:
            return getattr( self.data, self.var_enum_ )[:, self.cut_slider]
        if self.filament_slider_on == True:
            return  getattr( self.data, self.var_enum_ )[self.filament_slider, :]
        if self.cut_slider_on == False and self.filament_slider_on == False:
            return  getattr( self.data, self.var_enum_ )[ : , : ].flatten()

    traits_view = View( #VGroup( Item( 'data@', show_label=False ) ),
                        Item( 'var_enum', label='Variable' ),
                        HGroup( 
                        Item( 'cut_slider_on', label='Cut slider' ),
                        Item( 'cut_slider', show_label=False, springy=True,
                              enabled_when='cut_slider_on == True',
                              editor=RangeEditor( low_name='zero',
                                                    high_name='n_cuts',
                                                    mode='slider',
                                                    auto_set=False,
                                                    enter_set=False,
                                                    ),
                               ),
                        ),
                        HGroup( 
                        Item( 'filament_slider_on', label='Filament slider' ),
                        Item( 'filament_slider', show_label=False, springy=True,
                              enabled_when='filament_slider_on == True',
                              editor=RangeEditor( low_name='zero',
                                                    high_name='n_filaments',
                                                    mode='slider',
                                                    auto_set=False,
                                                    enter_set=False,
                                                    ),
                                )
                        ) )



if __name__ == '__main__':
    yarn = YMBSlider()
    yarn.configure_traits()


