

from enthought.traits.api import \
    HasTraits, Directory, List, Int, Float, Any, \
    on_trait_change, File, Constant, Instance, Trait, \
    Array, Str, Property, cached_property, WeakRef, \
    Dict, Button
    
from enthought.util.home_directory import \
    get_home_directory

from enthought.traits.ui.api import \
    View, Item, DirectoryEditor, TabularEditor, HSplit, VGroup, \
    TableEditor
    
from enthought.traits.ui.table_column import \
    ObjectColumn
    
from enthought.traits.ui.menu import \
    OKButton, CancelButton
    
from enthought.traits.ui.tabular_adapter \
    import TabularAdapter
    
import os

import csv

from numpy import array

from scipy.io import read_array
from numpy import loadtxt, argmax, polyfit, poly1d, frompyfunc

from enthought.traits.ui.table_filter \
    import EvalFilterTemplate, MenuFilterTemplate, RuleFilterTemplate, \
           EvalTableFilter


from enthought.traits.ui.file_dialog  \
    import open_file, FileInfo, TextInfo, ImageInfo


def fcomma2fdot(x): return float( replace(x, ',', '.') )

class InfoShellSpec(HasTraits):
    ''' Specification of the infograph shell-result set.
    
    Defines the parameters varied in the experiment design.
    The parameters can be read from the file included in the design.
    '''

    data_converters = { 0 : fcomma2fdot,
                        1 : fcomma2fdot }
                
    design_file    = Str( 'exdesign_num.csv')
    data_dir       = Str('data')
    data_file_name = " '%s %dnm %d.TRA' % (self.embedding, self.torque, self.rep_number ) "
    
class InfoShellReader(HasTraits):
    '''Read the data from the directory
    
    The design is described in semicolon-separated
    csv file providing the information about
    design parameters.
    
    Each file has the name n.txt
    '''

    #--------------------------------------------------------------------
    # Specification of the design - factor list, relative paths, etc
    #--------------------------------------------------------------------
    open_infoshell = Button()
    def _open_infoshell_fired(self):
        file_name = open_file( filter = ['*'], 
                               extensions = [FileInfo(), TextInfo()] )        
        if file_name != '':
            self.infoshell_spec_file = file_name
    
    infoshell_spec_file = File
    def _infoshell_spec_file_changed( self ):
        f = file( self.infoshell_spec_file )
        str = f.read()
        self.infoshell_spec = eval( 'InfoShellSpec( %s )' % str ) 

    infoshell_spec = Instance( InfoShellSpec )
    def _infoshell_spec_default(self):
        return InfoShellSpec()

    @on_trait_change('infoshell_spec')
    def _reset_design_file(self):
        dir = os.path.dirname( self. infoshell_spec_file )
        infoshell_file = self.infoshell_spec.design_file
        self.design_file =  os.path.join( dir, infoshell_file )

    #--------------------------------------------------------------------
    # file containing the association between the factor combinations 
    # and data files having the data 
    #--------------------------------------------------------------------
    design_file = File
    
    def _design_file_changed( self ):
        self.infoshell = self._read_infoshell()

    infoshell_table_columns = Property( List, depends_on = 'infoshell_spec+' )
    @cached_property
    def _get_infoshell_table_columns(self):
        return [ ObjectColumn( name = ps[2],
                               editable = False,
                               width = 0.15 ) for ps in self.infoshell_spec.factors ]

    infoshell = List( Any )
    def _infoshell_default(self):
        return self._read_infoshell()

    def _read_infoshell(self):
        ''' Read the experiment design. 
        '''
        if exists( self.design_file ):
            reader = csv.reader( open( self.design_file, 'r' ), delimiter=';' )
    
            data_dir = os.path.join( os.path.dirname( self.design_file ), 
                                     self.infoshell_spec.data_dir )
            
            return [ ExRun( self, row, data_dir = data_dir ) for row in reader ]
        else:
            return []

    selected_exrun = Instance( ExRun )
    def _selected_exrun_default(self):
        if len( self.infoshell ) > 0:
            return self.infoshell[0]
        else:
            return None

    last_exrun = Instance( ExRun )
    
    selected_exruns = List( ExRun )
    
    #------------------------------------------------------------------
    # Array plotting
    #-------------------------------------------------------------------
    # List of arrays to be plotted
    data = Instance(AbstractPlotData)
    def _data_default(self):
        return ArrayPlotData( x = array([]), y = array([]) )

    @on_trait_change('selected_exruns')
    def _rest_last_exrun(self):
        if len( self.selected_exruns ) > 0:
            self.last_exrun = self.selected_exruns[-1]
    
    @on_trait_change('selected_exruns')
    def _reset_data(self):
        '''
        '''
        runs, xlabels, ylabels, ylabels_fitted = self._generate_data_labels()
        for name in self.plot.plots.keys():
            self.plot.delplot( name )

        for idx, exrun in enumerate( self.selected_exruns ):
            if not self.plot.datasources.has_key( xlabels[idx] ):
                self.plot.datasources[ xlabels[idx] ] = ArrayDataSource( exrun.xdata,
                                                                         sort_order = 'none' )
            if not self.plot.datasources.has_key( ylabels[idx] ):
                self.plot.datasources[ ylabels[idx] ] = ArrayDataSource( exrun.ydata,
                                                                         sort_order = 'none' )

            if not self.plot.datasources.has_key( ylabels_fitted[idx] ):
                self.plot.datasources[ ylabels_fitted[idx] ] = ArrayDataSource( exrun.polyfit,
                                                                         sort_order = 'none' )

        for run, xlabel, ylabel, ylabel_fitted in zip( runs, xlabels, ylabels, ylabels_fitted ):
            self.plot.plot( (xlabel,ylabel), color = 'brown' )
            self.plot.plot( (xlabel,ylabel_fitted), color = 'blue' )

    def _generate_data_labels(self):
        ''' Generate the labels consisting of the axis and run-number.
        '''
        return ( map( lambda e: e.std_num, self.selected_exruns ),
                 map( lambda e: 'x-%d' % e.std_num, self.selected_exruns ),
                 map( lambda e: 'y-%d' % e.std_num, self.selected_exruns ),
                 map( lambda e: 'y-%d-fitted' % e.std_num, self.selected_exruns ) )
            
    
    plot = Instance(Plot)
    def _plot_default(self):
        p = Plot()
        p.tools.append(PanTool(p))
        p.overlays.append(ZoomTool(p))
        return p

    view_traits = View( HSplit( VGroup( Item('open_infoshell', 
                                             style = 'simple' ),
                                        Item('infoshell', 
                                             editor = exrun_table_editor, 
                                             show_label = False, style = 'custom' )
                                        ),
                                VGroup( Item('last_exrun@', 
                                             show_label = False ),
                                        Item('plot',
                                             editor=ComponentEditor(), 
                                             show_label = False,
                                             resizable = True
                                             ),
                                          ),
                                ),
#                        handler = InfoShellReaderHandler(),
                        resizable = True,
                        buttons = [ OKButton, CancelButton ],
                        height = 1.,
                        width = 1. )
    
    

doe_reader = InfoShellReader()
doe_reader.configure_traits( view = 'view_traits' )

