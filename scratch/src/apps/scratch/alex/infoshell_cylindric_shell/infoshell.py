from enthought.traits.api import \
    HasTraits, Directory, List, Int, Float, Any, \
    on_trait_change, File, Constant, Instance, Trait, \
    Array, Str, Property, cached_property, WeakRef, \
    Dict, Button
    
from enthought.util.home_directory import \
    get_home_directory

from enthought.traits.ui.api import \
    View, Item, DirectoryEditor, TabularEditor, HSplit, VGroup, \
    TableEditor, EnumEditor, Handler
    
from enthought.traits.ui.table_column import \
    ObjectColumn
    
from enthought.traits.ui.menu import \
    OKButton, CancelButton, Menu, MenuBar, Separator, Action
    
from enthought.traits.ui.tabular_adapter \
    import TabularAdapter
    
from string import strip
import os

import csv

from numpy import array

from scipy.io import read_array
from numpy import loadtxt, argmax, polyfit, poly1d, frompyfunc

from enthought.traits.ui.table_filter \
    import EvalFilterTemplate, MenuFilterTemplate, RuleFilterTemplate, \
           EvalTableFilter

from infoshell_file import InfoShellFile

from enthought.traits.ui.file_dialog  \
    import open_file, FileInfo, TextInfo, ImageInfo

import sys

def fcomma2fdot(x): return float( replace(x, ',', '.') )

#-----------------------------------------------------------------------------------
# ExDesignReader
#-----------------------------------------------------------------------------------
from enthought.enable.component_editor import \
    ComponentEditor
from enthought.chaco.api import \
    Plot, AbstractPlotData, ArrayPlotData, \
    ArrayDataSource
from enthought.chaco.tools.api import \
    PanTool, ZoomTool

infoshell_list_editor =  TableEditor(
    #columns_name = 'exdesign_table_columns',
    columns      = [ ObjectColumn( name = 'data_file',
                                   editable = False,
                                   width = 150 ) ],
    #selection_mode = 'rows',
    selected = 'object.selected_data_file',
#    auto_add     = False,
#    configurable = True,
#    sortable     = True,
#    sort_model  = True,   
    auto_size   = False,
    filters     = [ EvalFilterTemplate, MenuFilterTemplate, RuleFilterTemplate ],
    search      = EvalTableFilter()   
)        

class MFnWTHandler( Handler ):
    def export_data(self, info):
       sys.exit(0)

    def exit_file(self, info):
        sys.exit(0)

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
        file_name = open_file( filter = ['*.is'], 
                               extensions = [FileInfo(), TextInfo()] )        
        if file_name != '':
            self.infoshell_spec_file = file_name
    
    infoshell_spec_file = File
    def _infoshell_spec_file_changed( self ):
        f = file( self.infoshell_spec_file )
        
        self.infoshell_files = []
        while 1:
            line = f.readline()
            if not line:
                break
            data_file = InfoShellFile( data_file = strip( line ) )
            self.infoshell_files.append( data_file )
            
    infoshell_files = List
    # -----------------------------------------------------------------
    # GdG :
    # Print the selected data to a file with the same name as the 
    # raw input file with the ending '.outG'
    # -----------------------------------------------------------------
    print_to_file_gdg = Button()
    def _print_to_file_gdg_fired(self):
        # ------------------------------------------------------------
        # MX - Save the output-file:
        # ------------------------------------------------------------ 
        save_arr_MX_G = self.selected_data_file.MX_G
        save_list_MX_G = save_arr_MX_G.tolist()

        # specify the labels of the outputfile:
        label_list_MX_G = ['row_no', 'elem_no', 'node_no', \
                           'mx', 'my', 'mxy', 'nx' , 'ny', 'nxy', \
                           'm1', 'm2', 'alpha_M', \
                           'mx_M', 'my_M', 'nx_M', 'ny_M', \
                           'sig_n_MX', 'sig_m_MX', \
                           'eta_n_MX', 'eta_m_MX', 'eta_tot_MX'
                          ]
        
        # save the array input to file 'outG' applying formated strings
        data_file = self.selected_data_file.data_file
        output_GdG_file = open(data_file + '.outG','w')

        # print the input parameters defined
        output_GdG_file.write('------------------------- \n'+\
                              'Evaluated raw input file: \n'+\
                              '------------------------- \n '+\
                              '"' + str(data_file)+ '"' + '\n\n\n'\
                              '------------------------- \n'+\
                              'Defined input parameters: \n'+\
                              '------------------------- \n'+\
                              'thickness at border [m]: \n D_b [m] = %5.3f \n\n'\
                              'thickness at top [m]: \n D_t [m] = %5.3f \n\n'\
                              'tensile strength of the concrete[MPa]: \n f_ctk = %5.3f \n\n'\
                              'flexural tensile trength of the concrete [MPa]: \n f_m = %5.3f \n\n\n\n'\
                               %(self.selected_data_file.D_b, \
                                 self.selected_data_file.D_t, \
                                 self.selected_data_file.f_ctk, \
                                 self.selected_data_file.f_m
                                 ))
        
        #print results for GdG:
        caption_str_MX_G = str('GdG : Evaluation for case MX:')
        label_str_MX_G = str()
        dash_line_str = '\n'
        for i in range(len(label_list_MX_G)):
            label_str_MX_G = label_str_MX_G + '%010s '
            dash_line_str = dash_line_str + '---------- '
            
        header_str_MX_G = caption_str_MX_G + dash_line_str + '\n' + label_str_MX_G +  dash_line_str + '\n'
        
        output_GdG_file.write(header_str_MX_G %tuple(label_list_MX_G))
        
        for line in save_list_MX_G:
            output_GdG_file.write('%10.0f %10.0f %10.0f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f \n' \
                             %(line[0], line[1], line[2], line[3], line[4], line[5], line[6], line[7], line[8], line[9], line[10], line[11], line[12], line[13], line[14], line[15], line[16], line[17], line[18], line[19], line[20] ) )

        print 'Selected data has been printed to file "'+ data_file + '.outG"'
        output_GdG_file.close()

        # ------------------------------------------------------------
        # MY - Save the output-file:
        # ------------------------------------------------------------ 
        save_arr_MY_G = self.selected_data_file.MY_G
        save_list_MY_G = save_arr_MY_G.tolist()
        
        # specify the labels of the outputfile:
        label_list_MY_G = ['row_no', 'elem_no', 'node_no', \
                      'mx', 'my', 'mxy', 'nx' , 'ny', 'nxy', \
                      'm1', 'm2', 'alpha_M', 'mx_M', 'my_M', 'nx_M', 'ny_M', \
                      'sig_n_MY', 'sig_m_MY', \
                      'eta_n_MY', 'eta_m_MY', 'eta_tot_MY'
                      ]
        
        # save the array input to file 'outG' applying formated strings
        data_file = self.selected_data_file.data_file
        output_GdG_file = open(data_file + '.outG','a')
        
        caption_str_MY_G = str('GdG : Evaluation for case MY:')
        label_str_MY_G = str()
        dash_line_str = '\n'
        for i in range(len(label_list_MY_G)):
            label_str_MY_G = label_str_MY_G + '%010s '
            dash_line_str = dash_line_str + '---------- '
            
        header_str_MY_G = '\n \n \n' + caption_str_MY_G + dash_line_str + '\n' + label_str_MY_G +  dash_line_str + '\n'
        
        output_GdG_file.write(header_str_MY_G %tuple(label_list_MY_G))
        
        for line in save_list_MY_G:
            output_GdG_file.write('%10.0f %10.0f %10.0f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f \n' \
                             %(line[0], line[1], line[2], line[3], line[4], line[5], line[6], line[7], line[8], line[9], line[10], line[11], line[12], line[13], line[14], line[15], line[16], line[17], line[18], line[19], line[20] ) )
        
        output_GdG_file.close()
        # ------------------------------------------------------------
        # NX - Save the output-file:
        # ------------------------------------------------------------ 
        save_arr_NX_G = self.selected_data_file.NX_G
        save_list_NX_G = save_arr_NX_G.tolist()
        
        # specify the labels of the outputfile:
        label_list_NX_G = ['row_no', 'elem_no', \
                           'node_no', 'mx', \
                           'my', 'mxy', 'nx' , 'ny', 'nxy', \
                           'n1', 'n2', 'alpha_N', \
                           'mx_N', 'my_N', 'nx_N', 'ny_N', \
                           'sig_n_NX', 'sig_m_NX', \
                           'eta_n_NX', 'eta_m_NX', 'eta_tot_NX'
                           ]
        
        # save the array input to file 'outG' applying formated strings
        data_file = self.selected_data_file.data_file
        output_GdG_file = open(data_file + '.outG','a')
        
        caption_str_NX_G = str('GdG : Evaluation for case NX:')
        label_str_NX_G = str()
        dash_line_str = '\n'
        for i in range(len(label_list_NX_G)):
            label_str_NX_G = label_str_NX_G + '%010s '
            dash_line_str = dash_line_str + '---------- '
            
        header_str_NX_G = '\n \n \n' + caption_str_NX_G + dash_line_str + '\n' + label_str_NX_G +  dash_line_str + '\n'
        
        output_GdG_file.write(header_str_NX_G %tuple(label_list_NX_G))
        
        for line in save_list_NX_G:
            output_GdG_file.write('%10.0f %10.0f %10.0f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f \n' \
                             %(line[0], line[1], line[2], line[3], line[4], line[5], line[6], line[7], line[8], line[9], line[10], line[11], line[12], line[13], line[14], line[15], line[16], line[17], line[18], line[19], line[20] ) )
        
        output_GdG_file.close()
        # ------------------------------------------------------------
        # NY - Save the output-file:
        # ------------------------------------------------------------ 
        save_arr_NY_G = self.selected_data_file.NY_G
        save_list_NY_G = save_arr_NY_G.tolist()
        
        # specify the labels of the outputfile:
        label_list_NY_G = ['row_no', 'elem_no', 'node_no', \
                      'mx', 'my', 'mxy', 'nx', 'ny', 'nxy', \
                      'n1', 'n2', 'alpha_N', \
                      'mx_N', 'my_N', 'nx_N', 'ny_N', \
                      'sig_n_NY', 'sig_m_NY', \
                      'eta_n_NY', 'eta_m_NY', 'eta_tot_NY'
                      ]
        
        # save the array input to file 'outG' applying formated strings
        data_file = self.selected_data_file.data_file
        output_GdG_file = open(data_file + '.outG','a')
        
        caption_str_NY_G = str('GdG : Evaluation for case NY:')
        label_str_NY_G = str()
        dash_line_str = '\n'
        for i in range(len(label_list_NY_G)):
            label_str_NY_G = label_str_NY_G + '%010s '
            dash_line_str = dash_line_str + '---------- '
            
        header_str_NY_G = '\n \n \n' + caption_str_NY_G + dash_line_str + '\n' + label_str_NY_G +  dash_line_str + '\n'
        
        output_GdG_file.write(header_str_NY_G %tuple(label_list_NY_G))
        
        for line in save_list_NY_G:
            output_GdG_file.write('%10.0f %10.0f %10.0f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f \n' \
                             %(line[0], line[1], line[2], line[3], line[4], line[5], line[6], line[7], line[8], line[9], line[10], line[11], line[12], line[13], line[14], line[15], line[16], line[17], line[18], line[19], line[20] ) )
        
        output_GdG_file.close()
    # -----------------------------------------------------------------
    # GdT :
    # Print the selected data to a file with the same name as the 
    # raw input file with the ending '.outT'
    # -----------------------------------------------------------------
    print_to_file_gdt = Button()
    def _print_to_file_gdt_fired(self):
        # ------------------------------------------------------------
        # MX - Save the output-file:
        # ------------------------------------------------------------ 
        save_arr_MX_T = self.selected_data_file.MX_T
        save_list_MX_T = save_arr_MX_T.tolist()

        # specify the labels of the outputfile:
        label_list_MX_T = ['row_no', 'elem_no', 'node_no',\
                           'mx', 'my', 'mxy', 'nx' , 'ny', 'nxy', \
                           'm1', 'm2', 'alpha_M', \
                           'mx_M', 'my_M', 'nx_M', 'ny_M', \
                           'e_MX', 'm_Eds_MX', \
                           'f_t_MX', 'f_Rtex_MX', \
                           'n_MX'
                           ]
        
        # save the array input to file 'outG' applying formated strings
        data_file = self.selected_data_file.data_file
        output_GdT_file = open(data_file + '.outT','w')

        # print the input parameters defined
        output_GdT_file.write('------------------------- \n'+\
                              'Evaluated raw input file: \n'+\
                              '------------------------- \n '+\
                              '"' + str(data_file)+ '"' + '\n\n\n'\
                              '------------------------- \n'+\
                              'Defined input parameters: \n'+\
                              '------------------------- \n'+\
                              'thickness at border [m]: \n D_b [m] = %5.3f \n\n'\
                              'thickness at top [m]: \n D_t [m] = %5.3f \n\n'\
                              'tensile strength of the reinforcement (longitudinal) [kN/m]: \n F_Rtex_l = %5.3f \n'\
                              'tensile strength of the reinforcement (transversal) [kN/m]: \n F_Rtex_q = %5.3f \n\n\n\n'\
                               %(self.selected_data_file.D_b, \
                                 self.selected_data_file.D_t, \
                                 self.selected_data_file.F_Rtex_l, \
                                 self.selected_data_file.F_Rtex_q
                                 ))
        #print results for GdG:
        caption_str_MX_T = str('GdT : Evaluation for case MX:')
        label_str_MX_T = str()
        dash_line_str = '\n'
        for i in range(len(label_list_MX_T)):
            label_str_MX_T = label_str_MX_T + '%010s '
            dash_line_str = dash_line_str + '---------- '
            
        header_str_MX_T = caption_str_MX_T + dash_line_str + '\n' + label_str_MX_T +  dash_line_str + '\n'
        
        output_GdT_file.write(header_str_MX_T %tuple(label_list_MX_T))
        
        for line in save_list_MX_T:
            output_GdT_file.write('%10.0f %10.0f %10.0f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f \n' \
                             %(line[0], line[1], line[2], line[3], line[4], line[5], line[6], line[7], line[8], line[9], line[10], line[11], line[12], line[13], line[14], line[15], line[16], line[17], line[18], line[19], line[20] ) )

        output_GdT_file.close()
        # ------------------------------------------------------------
        # MY - Save the output-file:
        # ------------------------------------------------------------ 
        save_arr_MY_T = self.selected_data_file.MY_T
        save_list_MY_T = save_arr_MY_T.tolist()

        # specify the labels of the outputfile:
        label_list_MY_T = ['row_no', 'elem_no', 'node_no',\
                           'mx', 'my', 'mxy', 'nx' , 'ny', 'nxy', \
                           'm1', 'm2', 'alpha_M', \
                           'mx_M', 'my_M', 'nx_M', 'ny_M', \
                           'e_MY', 'm_Eds_MY', \
                           'f_t_MY', 'f_Rtex_MY', \
                           'n_MY'
                           ]
        
        # save the array input to file 'outG' applying formated strings
        data_file = self.selected_data_file.data_file
        output_GdT_file = open(data_file + '.outT','a')

        #print results for GdG:
        caption_str_MY_T = str('GdT : Evaluation for case MY:')
        label_str_MY_T = str()
        dash_line_str = '\n'
        for i in range(len(label_list_MY_T)):
            label_str_MY_T = label_str_MY_T + '%010s '
            dash_line_str = dash_line_str + '---------- '
            
        header_str_MY_T = '\n \n \n' + caption_str_MY_T + dash_line_str + '\n' + label_str_MY_T +  dash_line_str + '\n'
        
        output_GdT_file.write(header_str_MY_T %tuple(label_list_MY_T))
        
        for line in save_list_MY_T:
            output_GdT_file.write('%10.0f %10.0f %10.0f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f \n' \
                             %(line[0], line[1], line[2], line[3], line[4], line[5], line[6], line[7], line[8], line[9], line[10], line[11], line[12], line[13], line[14], line[15], line[16], line[17], line[18], line[19], line[20] ) )

        output_GdT_file.close()

        # ------------------------------------------------------------
        # NX - Save the output-file:
        # ------------------------------------------------------------ 
        save_arr_NX_T = self.selected_data_file.NX_T
        save_list_NX_T = save_arr_NX_T.tolist()

        # specify the labels of the outputfile:
        label_list_NX_T = ['row_no', 'elem_no', 'node_no',\
                           'mx', 'my', 'mxy', 'nx' , 'ny', 'nxy', \
                           'n1', 'n2', 'alpha_N', \
                           'mx_N', 'my_N', 'nx_N', 'ny_N', \
                           'e_NX', 'm_Eds_NX', \
                           'f_t_NX', 'f_Rtex_NX', \
                           'n_NX'
                           ]
        
        # save the array input to file 'outG' applying formated strings
        data_file = self.selected_data_file.data_file
        output_GdT_file = open(data_file + '.outT','a')

        #print results for GdG:
        caption_str_NX_T = str('GdT : Evaluation for case NX:')
        label_str_NX_T = str()
        dash_line_str = '\n'
        for i in range(len(label_list_NX_T)):
            label_str_NX_T = label_str_NX_T + '%010s '
            dash_line_str = dash_line_str + '---------- '
            
        header_str_NX_T = '\n \n \n' + caption_str_NX_T + dash_line_str + '\n' + label_str_NX_T +  dash_line_str + '\n'
        
        output_GdT_file.write(header_str_NX_T %tuple(label_list_NX_T))
        
        for line in save_list_NX_T:
            output_GdT_file.write('%10.0f %10.0f %10.0f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f\n' \
                             %(line[0], line[1], line[2], line[3], line[4], line[5], line[6], line[7], line[8], line[9], line[10], line[11], line[12], line[13], line[14], line[15], line[16], line[17], line[18], line[19], line[20] ) )

        output_GdT_file.close()

        # ------------------------------------------------------------
        # NY - Save the output-file:
        # ------------------------------------------------------------ 
        save_arr_NY_T = self.selected_data_file.NY_T
        save_list_NY_T = save_arr_NY_T.tolist()

        # specify the labels of the outputfile:
        label_list_NY_T = ['row_no', 'elem_no', 'node_no',\
                           'mx', 'my', 'mxy', 'nx' , 'ny', 'nxy', \
                           'n1', 'n2', 'alpha_N', \
                           'mx_N', 'my_N', 'nx_N', 'ny_N', \
                           'e_NY', 'm_Eds_NY', \
                           'f_t_NY', 'f_Rtex_NY', \
                           'n_NY'
                           ]
        
        # save the array input to file 'outG' applying formated strings
        data_file = self.selected_data_file.data_file
        output_GdT_file = open(data_file + '.outT','a')

        #print results for GdG:
        caption_str_NY_T = str('GdT : Evaluation for case NY:')
        label_str_NY_T = str()
        dash_line_str = '\n'
        for i in range(len(label_list_NY_T)):
            label_str_NY_T = label_str_NY_T + '%010s '
            dash_line_str = dash_line_str + '---------- '
            
        header_str_NY_T = '\n \n \n' + caption_str_NY_T + dash_line_str + '\n' + label_str_NY_T +  dash_line_str + '\n'
        
        output_GdT_file.write(header_str_NY_T %tuple(label_list_NY_T))
        
        for line in save_list_NY_T:
            output_GdT_file.write('%10.0f %10.0f %10.0f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f \n' \
                             %(line[0], line[1], line[2], line[3], line[4], line[5], line[6], line[7], line[8], line[9], line[10], line[11], line[12], line[13], line[14], line[15], line[16], line[17], line[18], line[19], line[20] ) )

        print 'Selected data has been printed to file "'+ data_file + '.outT"'
        output_GdT_file.close()

    selected_data_file = Instance( InfoShellFile )
    def _selected_data_file_default(self):
        if len( self.infoshell_files ) > 0:
            return self.infoshell_files[0]
        else:
            return None

    view_traits = View( HSplit( VGroup( Item('open_infoshell', 
                                             style = 'simple', show_label = False ),
                                        Item('print_to_file_gdg', label = 'Print to file (.outG)', style='simple', show_label = False),
                                        Item('print_to_file_gdt', label = 'Print to file (.outT)', style='simple', show_label = False),
                                        Item('infoshell_files',
                                             editor = infoshell_list_editor, show_label = False, width=100 ),
                                        ),
                                        Item( 'selected_data_file',
                                              show_label = False, style = 'custom',
                                              resizable = True,
                                              width = 800
                                        ),
                                ),
#                        handler = InfoShellReaderHandler(),
                        resizable = True,
                        buttons = [ OKButton, CancelButton ],
                        menubar=MenuBar(Menu(Action(name="Export (now works only as 'Exit')", action="export_data"),
                                             Action(name="E&xit", action="exit_file"),
                                             name = 'File')),
                        handler = MFnWTHandler,
                        height = 0.6,
                        width = 0.9 )
    
doe_reader = InfoShellReader( )
doe_reader.configure_traits( view = 'view_traits' )
