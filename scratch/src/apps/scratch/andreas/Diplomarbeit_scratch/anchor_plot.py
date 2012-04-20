from numpy import *
from matplotlib import pyplot

from enthought.traits.api import \
    HasTraits, Directory, List, Int, Float, Any, Enum, \
    on_trait_change, File, Constant, Instance, Trait, \
    Array, Str, Property, cached_property, WeakRef, \
    Dict, Button, Color, Bool, DelegatesTo, Callable

import csv
import os

from promod.simdb import \
    SimDB

from enthought.util.home_directory import \
    get_home_directory


simdb = SimDB()


class anchor_ouput( HasTraits ):

        def read_data( self, file_name ):
            '''read state data of the test for anchors, plates etc.
            '''
            print '*** read test data from file: %s ***' % ( file_name )

            file = open( file_name, 'r' )
            lines = file.readlines()

            for i in range( 0, len( lines ) ):
                for entries in lines[i].split( ';' ):
                    if  entries == '"Zeit sec"':
                        row_data = i + 1

            column_headings_arr = array( lines[row_data - 1].split( ';' ) )
            idx_u = where( column_headings_arr == '"Verformung mm"' )[0]
            idx_f = where( column_headings_arr == '"Kraft N"' )[0]

            idx_verformung = where
            input_arr = loadtxt( file_name , delimiter = ';', skiprows = row_data, dtype = str )
            shape = input_arr.shape
            for i in range ( 0, shape[0] ):
                for j in range ( 0, shape[1] ):
                    input_arr[i, j] = float( input_arr[i, j].replace( ",", "." ) )
            input_arr = array( input_arr, dtype = float )

            return {'f':abs( input_arr[:, idx_f] ) / 1000, 'u': abs( input_arr[:, idx_u] )}

        def read_data_spring( self, file_name ):
            '''read state data of the test for anchors, plates with Wegaufnehemer etc.
            '''
            print '*** read test data from file: %s ***' % ( file_name )

            file = open( file_name, 'r' )
            lines = file.readlines()

            input_arr = loadtxt( file_name , delimiter = ';', dtype = float )
            shift_idx = 0
            u = abs( input_arr[shift_idx:, 2] + input_arr[shift_idx:, 3] ) / 2 - abs( input_arr[shift_idx, 2] + input_arr[shift_idx, 3] ) / 2
            f = input_arr[shift_idx:, 1]
            return {'f': f , 'u': u}


        def set_plot( self ):
            self.fig = pyplot.figure( facecolor = "white" )
            self.ax1 = self.fig.add_subplot( 1, 1, 1 )
            self.ax1.set_xlabel( 'Verformung [mm]', fontsize = 22 )
            self.ax1.set_ylabel( 'Kraft [kN]', fontsize = 22 )
            self.ax1.legend

        def plot_file( self, x, y, filename = 'name' ):
            max_y = argmax( y ) + 20
            x = x[:max_y]
            y = y[:max_y]
            self.ax1.plot( x, y, label = str( filename ), linewidth = 1.5 )

        def save( self, savename ):
            self.ax1.legend( loc = 5 )
            self.ax1.grid( True )
            for xlabel_i in self.ax1.get_xticklabels():
                xlabel_i.set_fontsize( 15 )
            for ylabel_i in self.ax1.get_yticklabels():
                ylabel_i.set_fontsize( 15 )
            self.fig.savefig( savename, orientation = 'portrait', bbox_inches = 'tight' )
            pyplot.clf()


        def show( self ):
            self.ax1.legend( loc = 4 )
            pyplot.show()

if __name__ == '__main__':

    cls = anchor_ouput()
    data_dir = os.path.join( simdb.simdb_dir,
                            'exdata',
                            'anchor_pull_out_tests'
#                            'anchor_shear_tests'
                             )

    serie_name_list = [
#                 'PO-6g6c-gb6cm',

#                 'PO-6g6c-s6cm',

#                  'PO-12g-tp',

#                  'PO-16u-gb',

#                  'PO-16u-gb',

                  'PO-16u-a',


#                  'PO-16u-p',

#                  'ST-6g6c-p',

#                  'ST-16u-a',


#                  'ST-12g-s6cm',

#                  'ST-12g-tp',

#                  'ST-16g-gb6cm',

#                  'ST-16u-gb',
                  ]
    serie_name = serie_name_list [0]

    file_name_list = [
#                      'PO-6g6c-gb6cm-V1a',
#                      'PO-6g6c-gb6cm-V1b',

#                      'PO-6g6c-s6cm-V1a',
#                      'PO-6g6c-s6cm-V1b',

#                      'PO-12g-tp'

#                      'PO-16u-gb-V1',
#                      'PO-16u-gb-V2b',
#                      'PO-16u-gb-V3',

#                      'PO-16u-p-V1',
#                      'PO-16u-p-V2a',


#                      'ST-16u-a-V1',

#                      'PO-16u-gb-V3',

#                       'ST-6g6c-p-V1',
#                       'ST-6g6c-p-V2',

#                       'ST-12g-s6cm',

#                       'ST-12g-tp',

                      'PO-16u-a-V2',


#                        'ST-16g-gb6cm-V1a',
#                        'ST-16g-gb6cm-V1b',

#                        'ST-16u-gb-V1b',
#                        'ST-16u-gb-V2c',
                      ]


    legend_list = [
#                   'PO-6g6C-gb6cm',

#                    'PO-6g6c-s6cm-V1',

#                    'PO-12g-tp'

#                      'PO-16g-gb-V1',
#                      'PO-16g-gb-V2',
#                      'PO-16g-gb-V3',

#                      'PO-16g-p-V1',
#                      'PO-16g-p-V2',

#                       'ST-6g6C-p-V1',
#                       'ST-6g6C-p-V2',

#                     'ST-12g-s6cm',

#                     'ST-12g-tp',

#                      'ST-16g-gb6cm-V1a',
#                      'ST-16g-gb6cm',

#                        'ST-16g-gb-V1',
#                        'ST-16g-gb-V2',

                      'PO-16g-a',
#                      'ST-6g6c-p-V2' ,
#                      'ST-16g-a',
#                      'St-16u-a-V3'
#                      'PO-16u-gb-V3'
                      ]

#    file_name_list = [
#                      'PO-16u-gb-V1',
##                      'PO-16u-gb-V2a' ,
##                      'PO-16u-gb-V2b',
#                      'PO-16u-gb-V3']
#
#    legend_list = [
#                   'PO-16u-gb-V1',
##                   'PO-16u-gb-V2a',
##                   'PO-16u-gb-V2b',
#                   'PO-16u-gb-V2']

    cls.set_plot()
    for file_name, legend in zip( file_name_list, legend_list ):
        data_dict = cls.read_data( os.path.join( data_dir, serie_name, file_name ) + '.raw' )
#        data_dict = cls.read_data_spring( os.path.join( data_dir, serie_name, file_name ) + '.csv' )

        cls.plot_file( data_dict['u'], data_dict['f'], filename = legend )

    cls.save( serie_name, )
#    cls.show()

