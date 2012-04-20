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
    ones, transpose, shape, append, argmin, fabs, identity, \
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


class csv( HasTraits ):



    file_name_s = os.path.join( '/home', 'andreas', 'Documents', 'materialdaten', 'Schwinden.csv' )
    file_name_k = os.path.join( '/home', 'andreas', 'Documents', 'materialdaten', 'Kriechen.csv' )

    def read_data_s( self ):
        '''read state data and geo date from csv-file using ';' as filed delimiter and ' ' (blank)
        as text delimiter.
        '''
        file_name_s = self.file_name_s
        print '*** read shrinkage data from file: %s ***' % ( file_name_s )

        # get the column headings defined in the first row 
        # of the csv hinge force input file
        # column_headings = array(["elem_no", "N_ip", "V_ip", "V_op"])
        #
        file = open( file_name_s, 'r' )
        lines = file.readlines()
        column_headings = lines[0].split( ';' )
        file.close()

        input_array = loadtxt( file_name_s , delimiter = ';', skiprows = 2 )
        alter = input_array[:, 0]
        eps_s = input_array[:, 1]


        # plot
        #
        from matplotlib import rc

#        from matplotlib import rc
#        rc( 'text', usetex = True )
#        rc( 'font', **{'family':'serif', 'serif':'Times'} )
#        rc( 'text.latex',  )

        from matplotlib import font_manager
        prop = font_manager.FontProperties( size = 20 )

        x_label = 'Alter [d]'
        y_label = 'Dehnung  [$^o/_{oo}$]'


        fig = pyplot.figure( facecolor = "white" )
        ax1 = fig.add_subplot( 1, 1, 1 )
        ax1.plot( alter, eps_s, color = 'black', label = '$Schwinden$' , linewidth = 1.5 )

        x_ticks = list( alter[5:] )
        ax1.set_xticks( tuple( x_ticks ) )

        x_ticks_label = join( ["%s" % el for el in x_ticks] )
        ax1.set_xticklabels( x_ticks_label )

        for xlabel_i in ax1.get_xticklabels():
            xlabel_i.set_fontsize( 20 )

        ax1.set_xlabel( x_label, fontsize = 20 )
        ax1.set_ylabel( y_label, fontsize = 20 )

        ax1.legend( ['HEYX'] )
        pyplot.show()


    def read_data_k( self ):
        '''read state data and geo date from csv-file using ';' as filed delimiter and ' ' (blank)
        as text delimiter.
        '''
        file_name_k = self.file_name_k
        print '*** read creep data from file: %s ***' % ( file_name_k )

        # get the column headings defined in the first row 
        # of the csv hinge force input file
        # column_headings = array(["elem_no", "N_ip", "V_ip", "V_op"])
        #
        file = open( file_name_k, 'r' )
        lines = file.readlines()
        column_headings = lines[0].split( ';' )
        file.close()

        input_array = loadtxt( file_name_k , delimiter = ';', skiprows = 2 )
        alter = input_array[:, 0]
        delta_all = input_array[:, 1]
        delta_schwinden = input_array[:, 2]
        delta_elastisch = input_array[:, 3]
        delta_kriechen = input_array[:, 4]
        endkriechzahl = input_array[:, 4] / input_array[:, 3]




        # plot
        #
#        from matplotlib import rc
#        rc( 'text', usetex = True )
#        rc( 'font', **{'family':'serif', 'serif':'Times'} )
#        rc( 'text.latex',  )

        from matplotlib import font_manager
        prop = font_manager.FontProperties( size = 20 )

        x_label = 'Alter [d]'
        y_label = 'Dehnung $[^o/_{oo}]$'


        fig = pyplot.figure( facecolor = "white" )
        ax1 = fig.add_subplot( 1, 1, 1 )
#        ax1.plot( alter, delta_schwinden, color = 'black', label='schwinden', linewidth = 1.5 )
        ax1.plot( alter, delta_elastisch, 'k-.', label = 'elastsich', linewidth = 1.5 )
        ax1.plot( alter, delta_kriechen, color = 'black', label = 'kriechen', linewidth = 1.5 )
#        ax1.plot( alter, delta_kriechen + delta_kriechen +delta_schwinden , color= 'black', label='gesamt', linewidth = 1.5 )


        x_ticks = list( alter[6:] )
        ax1.set_xticks( tuple( x_ticks ) )
        ax1.set_ylim( ( -2.0, .0 ) )
        x_ticks_label = join( ["%s" % el for el in x_ticks] )
        ax1.set_xticklabels( x_ticks_label )
#
#        
#

        for tick in ax1.xaxis.get_major_ticks():
            tick.label1.set_fontsize( 15 )
        for tick in ax1.yaxis.get_major_ticks():
            tick.label1.set_fontsize( 15 )
        ax1.set_xlabel( x_label, fontsize = 22 )
        ax1.set_ylabel( y_label, fontsize = 22 )

        ax1.legend( prop = prop )
        pyplot.show()




if __name__ == '__main__':

    csv = csv()
    csv.read_data_s()
#    csv.read_data_k()
