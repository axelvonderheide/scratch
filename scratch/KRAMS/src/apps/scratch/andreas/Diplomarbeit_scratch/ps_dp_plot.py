'''
Created on Oct 19, 2010

@author: andreas
'''
from enthought.traits.api import \
    HasTraits, Directory, List, Int, Float, Any, Enum, \
    on_trait_change, File, Constant, Instance, Trait, \
    Array, Str, Property, cached_property, WeakRef, \
    Dict, Button, Color, Bool, DelegatesTo, Callable

from enthought.util.home_directory import \
    get_home_directory

from enthought.traits.ui.api import \
    View, Item, DirectoryEditor, TabularEditor, HSplit, Tabbed, VGroup, \
    TableEditor, Group, ListEditor, VSplit, HSplit, VGroup, HGroup, Spring, \
    Include

from enthought.traits.ui.table_column import \
    ObjectColumn

from enthought.traits.ui.menu import \
    OKButton, CancelButton

from enthought.traits.ui.tabular_adapter \
    import TabularAdapter

from numpy import \
    array, loadtxt, arange, sqrt, zeros, arctan, sin, cos, ones_like, \
    vstack, savetxt, hstack, argsort, fromstring, zeros_like, shape, \
    copy, c_, newaxis, argmax, where, sqrt, frompyfunc, sum, \
    ones, transpose, shape, append, argmin, fabs, identity, unique, \
    max as ndmax, min as ndmin


from math import pi
from string import split
import os

from scipy.io import read_array

from promod.simdb import \
    SimDB

import pickle
import string
from os.path import join

from matplotlib import pyplot


class ReadPS( HasTraits ):



    file_name = os.path.join( '/home/andreas/simdb/simdata/output_data_mushroof/displacement',
#                              'ps_verf_varA_w_asym_s_asym.csv',
                              'ps_verf_varA_g.csv'
                              )

    def read_data( self ):
        '''read state data and geo date from csv-file using ';' as filed delimiter and ' ' (blank)
        as text delimiter.
        '''
        file_name = self.file_name

        # get the column headings defined in the first row 
        # of the csv hinge force input file
        # column_headings = array(["elem_no", "N_ip", "V_ip", "V_op"])
        #
        file = open( file_name, 'r' )
        lines = file.readlines()
        column_headings = lines[0].split( ';' )
        file.close()
        from matplotlib import rc
        rc( 'text', usetex = True )
        rc( 'font', **{'family':'serif', 'serif':'times'} )
        rc( 'text.latex', )
        input_array = loadtxt( file_name , delimiter = ';', skiprows = 1 )

        param_1 = unique( input_array[:, 0] )
        param_2 = unique( input_array[:, 1] )

        fig = pyplot.figure( facecolor = "white" )
        ax1 = fig.add_subplot( 1, 1, 1 )


        for param in param_2:
            x = input_array[:, 0][where( input_array[:, 1] == param )]
            y = input_array[:, 2][where( input_array[:, 1] == param )]
            ax1.plot( x, y, label = str( param ) , linewidth = 1.5 )

        ax1.set_xlabel( '$w_{St,f}$[m]', fontsize = 20 )
        ax1.set_ylabel( '$u_{z}$[mm]', fontsize = 20 )

        ax1.grid()

        ax1.legend( loc = 0 )
        pyplot.show()







if __name__ == '__main__':

    ReadPS = ReadPS()
#    csv.read_data_s()
    ReadPS.read_data()
