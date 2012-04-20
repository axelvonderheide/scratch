'''
Created on 9.11.2010

@author: Vasek
'''
from enthought.traits.api import HasTraits, Float, Property, cached_property, \
    Event, Array, Instance, File, Int, Directory, Button, Range, Enum
from enthought.traits.ui.api import View, Item, FileEditor, DirectoryEditor, \
    HistoryEditor
from enthought.traits.ui.menu import OKButton, CancelButton, Action, Menu, \
    MenuBar
from mpl_toolkits.mplot3d import Axes3D, axes3d
from numpy import loadtxt, min, array, mean, std, arange, histogram, zeros_like, \
    meshgrid, ones_like, histogram2d, c_, equal, cumsum, vstack, hstack, savetxt, \
    sqrt, sum, all, zeros_like, zeros, ones
from os.path import join
from promod.simdb import SimDB
import matplotlib.pyplot as plt
import os
import re

simdb = SimDB()
data_dir = join( simdb.exdata_dir, 'trc', 'yarn_structure', 'VET', 'raw_data' )
#data_dir = join( simdb.exdata_dir, 'trc', 'yarn_structure', 'MAG', 'raw_data' )



class Data( HasTraits ):
    '''
        Read raw data into numpy arrays
    '''
    ddir = Directory( data_dir, auto_set=False, enter_set=True ) # data directory
    #n_cut = Int()
    
    slice = Property( depends_on='ddir' )
    @cached_property
    def _get_slice( self ):
        """
            Prepare data of continuous fibers 
        """
        slice_point_list, slice_radius_list, slice_len_list = self.read_slice()
        slice_len_arr = array( slice_len_list )
        print 'slice lens', slice_len_arr
        offset_arr = cumsum( slice_len_arr )
        slice_offset_arr = zeros_like( offset_arr )
        slice_offset_arr[1:] = offset_arr[:-1]
        print 'slice offsets', slice_offset_arr
    
        filam_connect_arr = self.connectivity
        # filter only continuous fibers
        filam_connect_arr = filam_connect_arr[ min( filam_connect_arr, axis=1 ) != -1. ]
        print filam_connect_arr.shape
    
        print filam_connect_arr.shape
        print slice_offset_arr.shape
    
        fil_map = array( filam_connect_arr + slice_offset_arr, dtype='int' )
    
        points = vstack( slice_point_list )
        radius = hstack( slice_radius_list )
    
        print points.shape
        print max( fil_map.flatten() )
    
        p = points[ fil_map.flatten() ]
        r = radius[ fil_map.flatten() ]

        return p, r
    
    slice_discontinuous = Property( depends_on='ddir' )
    @cached_property
    def _get_slice_discontinuous( self ):
        """
            Prepare data of continuous fibers 
        """
        slice_point_list, slice_radius_list, slice_len_list = self.read_slice()
        slice_len_arr = array( slice_len_list )
        print 'slice lens', slice_len_arr
        offset_arr = cumsum( slice_len_arr )
        slice_offset_arr = zeros_like( offset_arr )
        slice_offset_arr[1:] = offset_arr[:-1]
        print 'slice offsets', slice_offset_arr
    
        filam_connect_arr = self.connectivity
        # filter only continuous fibers
        filam_connect_arr = filam_connect_arr[ min( filam_connect_arr, axis=1 ) == -1. ]
        filam_connect_arr = filam_connect_arr[ all( filam_connect_arr, axis=1 ) != -1. ]
        print filam_connect_arr.shape
    
        print filam_connect_arr.shape
        print slice_offset_arr.shape
        fil_map = -ones_like( filam_connect_arr )
        slice_offset_arr = array( zeros_like( filam_connect_arr ) + slice_offset_arr )
        fil_map[filam_connect_arr > -1] = array( filam_connect_arr[filam_connect_arr > -1] + slice_offset_arr[filam_connect_arr > -1], dtype='int' ) 
        print 'fil map', fil_map
        fil_map = array( fil_map, dtype='int' )
        #fil_map = array( filam_connect_arr + slice_offset_arr, dtype='int' )
    
        points = vstack( [vstack( slice_point_list ), -ones( 3 )] ) 
        radius = hstack( [hstack( slice_radius_list ), -ones( 1 )] )
    
        print points.shape
        print max( fil_map.flatten() )
    
        p = points[ fil_map.flatten() ]
        r = radius[ fil_map.flatten() ]

        return p, r
    
    def read_slice( self ):
        """
            Given a directory name, returns a list of all slice files and slice numbers.
            Slice filenames must satisfy the following format '*Schnitt(num).txt'
            Return:
                x,y,z list of np arrays
                d list of np arrays
        """
        files = []
        num = []
        paths = os.listdir( self.ddir )  # list of paths in that dir
        for fname in paths:
            match = re.search( r'\w+Schnitt(\d+).txt', fname ) # find files of given pattern
            if match:
                files.append( os.path.abspath( os.path.join( self.ddir, fname ) ) )
                num.append( int( match.group( 1 ) ) )
        num, files = zip( *sorted( zip( num, files ) ) ) 
        slice_point_list = []
        slice_radius_list = []
        slice_len_list = []
        # TODO: load distance into slice distance array from slice_height
        slice_distance = [500] * len( num ) # micrometers
        for n, file, dist in zip( num, files, slice_distance ):
            print 'reading slice data_file -- slice %i' % n
            points = loadtxt( file ,
                              skiprows=1,
                              usecols=( 1, 2, 3 ) )
    
            y = points[ :, 0]
            z = points[ :, 1]
            x = ones_like( y ) * n * dist
            r = points[ :, 2]
    
            slice_point_list.append( c_[ x, y, z ] )
            slice_radius_list.append( r )
            slice_len_list.append( points.shape[0] )
        return slice_point_list, slice_radius_list, slice_len_list
    
    # TODO: get slice height
    slice_height = Property( depends_on='ddir' )
    @cached_property
    def _get_slice_height( self ):
        """
            Read data with cut height
            First cut has x-coordinate = 0
            Return np array
        """
        pass    

    connectivity = Property( depends_on='ddir' )
    @cached_property
    def _get_connectivity( self ):
        """
            Read connectivity data.
            Return np array
        """
        print 'reading data_file -- connectivity.txt'
        connectivity_file = join( self.ddir, 'connectivity.txt' )
        connect = loadtxt( connectivity_file, dtype='float', skiprows=1 )
        # replace 0 by -1 (not identified)
        connect[ equal( 0, connect )] = -1
        return connect
    
    bq = Property( Array, depends_on='ddir' )
    @cached_property
    def _get_bq( self ):
        """
            Read bond_quality data.
            Return np array
        """
        print 'reading data_file -- bond_quality.txt'
        bond_quality_file = join( self.ddir, 'bond_quality.txt' )
        bond_quality = loadtxt( bond_quality_file, dtype='float', skiprows=1, usecols=None )
        # filter continuous fibers
        # and the first number column
        return bond_quality[:, 1:][min( bond_quality[:, 1:], axis=1 ) != -1]
    
              
    view = View( Item( 'ddir',
                       id='ddir',
                       label='Directory',
                       editor=DirectoryEditor( entries=10 )
                ),
                )
    
    
    
    # read data from text files
    def __call__( self ):
        self._get_connectivity(), self._get_slice(), self._get_bq()
        



class Statistics( HasTraits ):
    
    #data = Data()
    
    mean = Property()
    def _get_mean( self, x ):
        """
            return mean value
            input numpy array 
        """
        return mean( x.T, axis=0 )
    
    # TODO: general, no bq
    def get_hist( self, x, bins=20 ):
            plt.hist( x, bins, normed=1 )
            

class Methods( HasTraits ):
    data = Instance( Data )
#    var_list = Enum( ['bq', 'bfl', 'slack', 'radius'] )
#    var_dict = {'bq' : data.bq,
#                'radius' : data.slice[1], }
    # TODO: settings of method
    method_list = Enum( ['yarn', 'slice'] )
    method_dict = {'yarn' :' data',
                   'slice' : 'get_slice', }

    def _get_length_between_cuts( self ):
        """
            Return np array of fiber length between cuts -- linear
        """
        p = self.data.slice[0]
        length = sqrt( sum( ( p[1:] - p[:-1] ) ** 2, axis=1 ) )
        print length
        return length
    
    def _get_slack( self ):
        return 0
    
    def _get_bond_free_length( self ):
        pass
    
    def get_bq_yarn( self ):
        return data.bq.flatten()
    
    #ix = Range()
    def get_bq_slice( self, ix ):
        '''
            return bond quality vector of the ix-th slice of the 
            experimental data
        '''
        return data.bq[:, ix]
    
        


if __name__ == '__main__':
#    stat = Statistics()
#    stat.get_bq_hist( 15 )
#    plt.show()
    data = Data()
    data.configure_traits()
    #data.set( ddir=data_dir )
    stats = Statistics( data=data )
    met = Methods( data=data )
    stats.get_hist( met.get_bq_yarn(), 15 )
    
    met._get_length_between_cuts()
    plt.show()
    print 'exit'
    print data.slice_discontinuous
    
