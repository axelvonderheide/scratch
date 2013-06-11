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

from enthought.traits.api import \
    HasTraits, Directory, List, Int, Float, Any, \
    on_trait_change, File, Constant, Instance, Trait, \
    Array, Str, Property, cached_property, WeakRef, \
    Dict, Button, Bool, Enum

from enthought.util.home_directory import \
    get_home_directory

from enthought.traits.ui.api import \
    View, Item, DirectoryEditor, TabularEditor, HSplit, VGroup, \
    TableEditor, EnumEditor, Handler
    
from enthought.traits.ui.table_column import \
    ObjectColumn
    
from enthought.traits.ui.menu import \
    OKButton, CancelButton
    
from enthought.traits.ui.tabular_adapter \
    import TabularAdapter
    
import os

import csv

from numpy import array, fabs, where, copy

from scipy.io import read_array
from numpy import loadtxt, argmax, polyfit, poly1d, frompyfunc

from enthought.traits.ui.table_filter \
    import EvalFilterTemplate, MenuFilterTemplate, RuleFilterTemplate, \
           EvalTableFilter

#-- Tabular Adapter Definition -------------------------------------------------

from string import replace
from os.path import exists

def fcomma2fdot(x): return float( replace(x, ',', '.') )

class ExDesignSpec(HasTraits):
    ''' Specification of the experiment design.
    
    Defines the parameters varied in the experiment design.
    The parameters can be read from the file included in the design.
    '''

    factors = []
    data_converters = { 0 : fcomma2fdot,
                        1 : fcomma2fdot }
                
    design_file    = Str('exdesign_num.csv')
    data_dir       = Str('data')
    data_file_name = " '%s %dnm %d.TRA' % (self.embedding, self.torque, self.rep_number ) "
    

#-----------------------------------------------------------------------------------
# ExDesignReader
#-----------------------------------------------------------------------------------
from enthought.enable.component_editor import \
    ComponentEditor
from enthought.chaco.api import \
    Plot, AbstractPlotData, ArrayPlotData, \
    ArrayDataSource, PlotAxis
from enthought.chaco.tools.api import \
    PanTool, ZoomTool

from enthought.traits.ui.file_dialog  \
    import open_file, FileInfo, TextInfo, ImageInfo

from enthought.traits.ui.api \
    import View, Item, TabularEditor
    
from enthought.traits.ui.tabular_adapter \
    import TabularAdapter
   
class ArrayAdapter ( TabularAdapter ):

    columns = Property
    def _get_columns(self):
        array_attr = getattr( self.object, self.name )
        if len( array_attr.shape ) == 1 :
            return []
        else:
            print 'resetting columns'
            n_columns = getattr( self.object, self.name ).shape[1]
            cols = [ (str(i), i ) for i in range( n_columns ) ]
            return [ ('i', 'index') ] + cols
        
    font        = 'Courier 10'
    alignment   = 'right'
    format      = '%6.2f'
    index_text  = Property
    
    def _get_index_text ( self ):
        return str( self.row )

array_editor = TabularEditor( adapter = ArrayAdapter() )
  
enum_dict = {}
 
class EXDesignReader(HasTraits):
    '''Read the data from the directory
    
    The design is described in semicolon-separated
    csv file providing the information about
    design parameters.
    
    Each file has the name n.txt
    '''
    #--------------------------------------------------------------------
    # Specification of the design - factor list, relative paths, etc
    #--------------------------------------------------------------------
    open_design = Button()
    def _open_design_fired(self):
        file_name = open_file( filter = ['*.DAT'], 
                               extensions = [FileInfo(), TextInfo()] )        
        if file_name != '':
            self.design_file = file_name
            self._reset_data()
    
    # File with description of the attributes
    exdesign_spec_file = File
    def _exdesign_spec_file_changed( self ):
        print 'changed file'
        f = file( self.exdesign_spec_file )
        str = f.read()
        self.exdesign_spec = eval( 'ExDesignSpec( %s )' % str ) 

    exdesign_spec = Instance( ExDesignSpec )
    def _exdesign_spec_default(self):
        print 'ExDesignSpec()', ExDesignSpec()
        return ExDesignSpec()

    @on_trait_change('exdesign_spec')
    def _reset_design_file(self):
        dir = os.path.dirname( self. exdesign_spec_file )
        exdesign_file = self.exdesign_spec.design_file
        self.design_file =  os.path.join( dir, exdesign_file )

    exdesign_table_columns = Property( List, depends_on = 'exdesign_spec+' )
    @cached_property
    def _get_exdesign_table_columns(self):
        return [ ObjectColumn( name = ps[2],
                               editable = False,
                               width = 0.15 ) for ps in self.exdesign_spec.factors ]

    #--------------------------------------------------------------------
    # file containing the association between the factor combinations 
    # and data files having the data 
    #--------------------------------------------------------------------
    design_file = File
    
    def _design_file_changed( self ):
        print 'reading design'
        self.exdesign = self._read_exdesign()

    exdesign = Array( float )
    def _exdesign_default(self):
        return self._read_exdesign()

    names_units_list = List
    exdesign_ironed = Array( float )
    exdesign_ironed_scaled = Array( float )
    jump_atol = Float

    def get_names_and_units(self, file_name):
        ''' Read the data from the source file.  '''
        names_units = []
        file = open(file_name,'r')    
        n = file.read().split()
        lenght = len(n)
        i = 0
        
        while i < lenght:
            if n[i] == '#BEGINCHANNELHEADER':
                desc = n[i+1]
                x = n[i+1].split(',')
                y = n[i+3].split(',')
                name_str = x[1] + ' [' + y[1] + ']' 
                names_units.append(name_str)
            i = i +1
        return names_units
    
    def _read_exdesign(self):
        ''' Read the experiment design. 
        '''
        if exists( self.design_file ):
            self.names_units_list = self.get_names_and_units(self.design_file)
            
            for i in range(len(self.names_units_list)):
                enum_dict[ self.names_units_list[i] ] = i
                
            self._reset_enum()
            
            ''' change the file name dat with asc.  '''  
            file_split = self.design_file.split('.')

            file_name = file_split[0] + '.ASC' 
            exdesign = loadtxt( file_name,  
                                delimiter = ';' )
            
            return exdesign
        else:
            return []
    
    @on_trait_change('jump_rtol, jump_atol')
    def _curve_ironing(self):
        '''removes the jumps in the curve due to 
        resetting of the displacement gauges '''
        
        print '*** curve ironing activated ***'
        
        # slice the column from the data array 
        # each column corresponds to a measured parameter 
        # e.g. displacement at a given point as function of time u = f(t))
        # NOTE: the first column must contain monotonic increasing time
        # this column is not ironed
        for idx in range(self.exdesign.shape[1]-1):

            # use copy in order to save initial state of the input data
            data_arr = copy(self.exdesign[:,idx+1])
            
            # get the difference between each point and its successor
            jump_arr =  data_arr[1:] - data_arr[0:-1]
            
            # get the value of the maximum relative jump
            jump_arr_max = max( fabs( jump_arr ))
            
            # get the maximum value of the measured data 
            data_arr_max = max( data_arr )
            data_arr_min = min( data_arr )

            data_arr_range = data_arr_max - data_arr_min
    
            # set the relative and absolute tolerances for the 
            # jump-criteria
            jump_rtol = self.jump_rtol 
            
            # determine the relevant criteria for a jump
            # based on the relative and absolute bound:
            jump_crit = max( self.jump_atol, jump_rtol * jump_arr_max, 
                             jump_rtol * data_arr_range ) 

            # get the indexes in 'data_column' after which a 
            # jump exceeds the defined tolerance criteria
            jump_indexes = where( fabs(jump_arr) > jump_crit )

            # glue the curve at each jump together
            for n in range(jump_indexes[0].shape[0]):
                # get the offsets at each jumb of the curve
                jidx = jump_indexes[0][n]
                shift = data_arr[jidx+1] - data_arr[jidx]

                # shift all succeeding values by the calculated offset
                data_arr[jidx+1:] -= shift
            print 'data_arr[:]= ', data_arr[:]
            self.exdesign_ironed[:,idx+1] = data_arr[:]
            
    @on_trait_change('scale_x, scale_y')
    def _curve_scaling(self):
        '''scale the values of the selected column with a scalar factor'''
        
        print '*** curve scaling activated ***'
        self.exdesign_ironed_scaled = copy(self.exdesign_ironed)
        
        x_index = float(self.xidx)
        y_index = float(self.yidx)

#        scale_x = float(self.scale_x)
#        scale_y = float(self.scale_y)

        print 'xxx1: self.exdesign_ironed_scaled[ :, x_index ]', self.exdesign_ironed_scaled[ :, x_index ]
        self.exdesign_ironed_scaled[ :, x_index ] = self.exdesign_ironed_scaled[ :, x_index ] * self.scale_x
        self.exdesign_ironed_scaled[ :, y_index ] = self.exdesign_ironed_scaled[ :, y_index ] * self.scale_y
        
        print 'xxx2: self.exdesign_ironed_scaled[ :, x_index ]', self.exdesign_ironed_scaled[ :, x_index ]
        print 'self.scale_x', self.scale_x
        print 'self.scale_y', self.scale_y
        
    plot = Instance(Plot)
    #plot = Instance(PlotAxis)
    def _plot_default(self):
        p = Plot()
        #p = PlotAxis()
        p.tools.append(PanTool(p))
        p.overlays.append(ZoomTool(p))
        return p

    xidx_enum = List( Str)    
    yidx_enum = List( Str)  

    xidx_name = Str 
    yidx_name = Str
    
    xidx = Str 
    yidx = Str


    def _xidx_name_changed(self): 
        print 'xidx changed'
        self.xidx = str( enum_dict.get(self.xidx_name) )
                    
    def _yidx_name_changed(self): 
        self.yidx = str( enum_dict.get(self.yidx_name) )
#        self._xidx_name_changed()
        
    jump_rtol = Float(0.03, auto_set = False, enter_set = True )
    scale_x = Float(1.0)   
    scale_y = Float(1.0)   

    def _reset_enum(self):
        self.xidx_enum = enum_dict.keys()
        self.yidx_enum = enum_dict.keys()
        self.xidx_name = enum_dict.keys()[0]
        self.yidx_name = enum_dict.keys()[0]

    @on_trait_change('xidx, yidx, scale_x, scale_y')
    def _reset_data(self):
        '''
        '''
        print '***RESET DATA***'
#       for scaling: look at auto_ticks, PlotAxis  
            
        for name in self.plot.plots.keys():
            self.plot.delplot( name )
        
        if self.exdesign_ironed.shape[0] == 0: 
            self.exdesign_ironed = copy(self.exdesign)
            self.exdesign_ironed_scaled = copy(self.exdesign)
        
        for idx in range( self.exdesign.shape[1] ):
            self.plot.datasources[ `idx` ] = ArrayDataSource( 
                                                self.exdesign_ironed_scaled[:,idx],
                                                sort_order = 'none')
        
        self.plot.plot( (self.xidx,self.yidx), color = 'black')
        self.plot.plot( (self.xidx,self.yidx), color = 'black')
        
        colors = ['red','blue','brown','yellow','black','orange','black',
                  'black','black','black','black']

    view_traits = View( HSplit( 
                            VGroup( Item('open_design',  
                                        show_label = False,
                                        style = 'simple' ),
                                    Item( name='xidx_name', 
                                          editor=EnumEditor(name='xidx_enum')),
                                    Item( name='yidx_name', 
                                          editor=EnumEditor(name='yidx_enum')),
                                    Item('jump_rtol'),
                                    Item('jump_atol'),
                                    Item('scale_x'),
                                    Item('scale_y'),
                                   ),
                                    Item('plot', 
                                        editor=ComponentEditor(), 
                                        show_label = False,
                                        resizable = True
                                        ),
                                ),
                        resizable = True,
                        buttons = [ OKButton, CancelButton ],
                        height = 0.8,
                        width = 0.8 )

doe_reader = EXDesignReader()
doe_reader.configure_traits( view = 'view_traits' )
