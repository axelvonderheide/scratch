'''
Created on Apr 30, 2010

@author: alexander
'''

#-------------------------------------------------------------------------------
#
# Copyright (c) 2009, IMB, RWTH Aachen.
# All rights reserved.
#
# This software is provided without warranty under the terms of the BSD
# license included in simvisage/LICENSE.txt and may be redistributed only
# under the conditions described in the aforementioned license.  The license
# is also available online at http://www.simvisage.com/licenses/BSD.txt
#
# Thanks for using Simvisage open source!
#
# Created on Aug 3, 2009 by: rch
#
# How to introduce ExperimentType class - as an class corresponding to SimModel
#
# - this class gathers the information about the inputs and outputs
#   of the object.
#
# - Each directory contains the `ex_type.cls' file giving the name of hte
#   ExperimentType subclass defining the inputs and outputs of the experiment
#   and further derived data needed for processing of experiments.
#
# - The ExperimentDB has a root directory and starts by scanning the subdirectories
#   for the `ex_type.cls' files. It verifies that the ExTypes are defined and known 
#   classes. Then, the mapping between a class and between the directories 
#   is established. Typically, there is a data pool in the home directory 
#   or a network file system accessible.
# 
# - The result of the scanning procedure is list of data directories available for 
#   each experiment type. Typically, data from a single treatment are grouped 
#   in the single directory, but this does not necessarily need to be the case. 
#   Therefore, the grouping is done independently on the data based on 
#   the values of the input factors.
#
#   Tasks: define the classes
#     - ETCompositeTensileTest 
#     - ETPlateTest
#     - ETPullOutTest, 
#     - ETYarnTensileTest
#
#   They inherit from class ExType specifying the inputs and outputs. They also specify
#   the default association of tracers to the x- and y- axis (plot templates). Further, grouping 
#   of equivalent values can be provided. The values of the inputs are stored 
#   in the directory in the ExType.db file using the pickle format. 
#
#   The ExTypeView class communicates with the particular ExType class. It should be able
#   to accommodate several curves within the plotting window. It should be possible to
#   select multiple runs to be plotted.
#      

# - Database of material components is stored in a separate directory tree. 
#   The identification
#   of the material components is provided using a universal key 
#   for a material component The component must be declared either as a matrix 
#   or reinforcement. It must be globally accessible from within 
#   the ExpTools and SimTools. 
#      

from enthought.traits.api import \
    HasTraits, Directory, List, Int, Float, Any, \
    on_trait_change, File, Constant, Instance, Trait, \
    Array, Str, Property, cached_property, WeakRef, \
    Dict, Button, Bool, Enum, Event, implements, DelegatesTo

from enthought.traits.ui.api import \
    View, Item, DirectoryEditor, TabularEditor, HSplit, VGroup, \
    TableEditor, EnumEditor, Handler, FileEditor, VSplit, Group

from enthought.traits.ui.menu import \
    OKButton, CancelButton

from enthought.util.home_directory import \
    get_home_directory

from enthought.traits.ui.table_column import \
    ObjectColumn

from enthought.traits.ui.tabular_adapter \
    import TabularAdapter

from enthought.traits.ui.table_filter import \
    EvalFilterTemplate, MenuFilterTemplate, RuleFilterTemplate, \
    EvalTableFilter

import os
import csv

from numpy import \
    array, fabs, where, copy, ones

from numpy import \
    loadtxt, argmax, polyfit, poly1d, frompyfunc, dot

#-- Tabular Adapter Definition -------------------------------------------------

from string import \
    replace

from os.path import \
    exists

import pickle

#-----------------------------------------------------------------------------------
# ExDesignReader
#-----------------------------------------------------------------------------------
from enthought.traits.ui.file_dialog  \
    import open_file, FileInfo, TextInfo, ImageInfo

class ExRun( HasTraits ):
    '''Read the data from the DAT file containing the measured data.
    
    The data is described in semicolon-separated
    csv file providing the information about
    data parameters.
    
    '''
    #--------------------------------------------------------------------
    # file containing the association between the factor combinations 
    # and data files having the data 
    #--------------------------------------------------------------------
    data_file = File

    # Derive the path to the file specifying the type of the experiment
    # The ex_type.cls file is stored in the same directory
    #
    ex_type_file_name = Property( depends_on = 'data_file' )
    @cached_property
    def _get_ex_type_file_name( self ):
        dir_path = os.path.dirname( self.data_file )
        file_name = os.path.join( dir_path, 'ex_type.cls' )
        return file_name

    # Derive the path to the pickle file storing the input data and derived output 
    # data associated with the experiment run 
    # 
    pickle_file_name = Property( depends_on = 'data_file' )
    @cached_property
    def _get_pickle_file_name( self ):
        return 'exrun.pickle'

    # Instance of the specialized experiment type with the particular
    # inputs and data derived outputs. 
    #
    ex_type = Instance( CCS )
    def _ex_type_default( self ):
        return CCS()

    def __init__( self, **kw ):
        '''Initialization: the ex_run is defined by the 
        data_file. The additional data - inputs and derived outputs
        are stored in the data_file.pickle. If this file exists, 
        the exrun is constructed from this data file.
        '''
        super( ExRun, self ).__init__( **kw )

        self.ex_type = CCS()
        self.unsaved = True

    def save_pickle( self ):
        '''Store the current state of the ex_run.
        '''
        file = open( self.pickle_file_name, 'w' )
        pickle.dump( self.ex_type, file )
        file.close()
        self.unsaved = False

    # Event to keep track of changes in the ex_type instance.
    # It is defined in order to inform the views about a change
    # in some input variable.
    #
    change_event = Event

    # Boolean variable set to true if the object has been changed. 
    # keep track of changes in the ex_type instance. This variable
    # indicates that the run is unsaved. The change may be later
    # canceled or confirmed.
    #
    unsaved = Bool( transient = True )

    @on_trait_change( 'ex_type.input_change' )
    def _set_changed_state( self ):
        print '***received input change event'''
        self.change_event = True
        self.unsaved = True


if __name__ == '__main__':
    exrun = ExRun()
    print 'pickle_file_name', exrun.pickle_file_name


    exrun.save_pickle()

    #exrun.configure_traits()

