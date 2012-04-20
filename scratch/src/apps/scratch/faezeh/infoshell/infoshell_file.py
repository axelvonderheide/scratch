
from enthought.traits.api import \
    HasTraits, Directory, List, Int, Float, Any, \
    on_trait_change, File, Constant, Instance, Trait, \
    Array, Str, Property, cached_property, WeakRef, \
    Dict, Button
    
from enthought.util.home_directory import \
    get_home_directory

from enthought.traits.ui.api import \
    View, Item, DirectoryEditor, TabularEditor, HSplit, Tabbed, VGroup, HSplit,\
    TableEditor, Group, ListEditor, VSplit
    
from enthought.traits.ui.table_column import \
    ObjectColumn
    
from enthought.traits.ui.menu import \
    OKButton, CancelButton
    
from enthought.traits.ui.tabular_adapter \
    import TabularAdapter
    
from numpy import array, loadtxt, arange, sqrt, zeros, arctan, sin, cos, ones_like,\
                vstack, savetxt, hstack, argsort, fromstring    

from math import pi
from string import split 
import os

#import csvPixel.py

from numpy import array, ones_like, arange

from scipy.io import read_array
from numpy import loadtxt, argmax, polyfit, poly1d, frompyfunc, sqrt

from enthought.traits.ui.table_filter \
    import EvalFilterTemplate, MenuFilterTemplate, RuleFilterTemplate, \
           EvalTableFilter

from enthought.pyface.dock.dock_sizer \
    import DockSizer, DockControl, DockRegion, DockStyle, DockSplitter, \
           no_dock_info, clear_window, features

#-- Tabular Adapter Definition -------------------------------------------------

class ArrayAdapterN ( TabularAdapter ):

    columns = [ ( 'i', 'index' ), ( 'Row.No.', 0 ), ( 'Elem.No.', 1 ),  ( 'Node.No.', 2 ), ('mx', 3), \
                ('my',4),         ('mxy',5),        ('nx',6),('ny',7),('nxy',8),('n1',9),('n2',10),('alpha_N',11),\
               ('mx_N',12),('my_N',13),('nx_N',14),('ny_N',15),  \
               ('sig_n',16),('sig_m',17),('eta_n',18),('eta_m',19), ( 'eta.tot', 16 ) ]
                
    font        = 'Courier 10'
    alignment   = 'right'
    format      = '%g'
   # column      = 'Elem.No.'
   # format      = '%d'
    index_text  = Property
#    index_image = Property
    
    def _get_index_text ( self ):
        return str( self.row )
        
    def _get_index_image ( self ):
        x, y, z = self.item
        if sqrt( (x - 0.5) ** 2 + (y - 0.5) ** 2 + (z - 0.5) ** 2 ) <= 0.25:
            return '@icons:red_ball'
            
        return None
 
class ArrayAdapterM ( TabularAdapter ):

    columns = [ ( 'i', 'index' ), ( 'Row.No.', 0 ), ( 'Elem.No.', 1 ),  ( 'Node.No.', 2 ), ('mx', 3), \
                ('my',4),         ('mxy',5),        ('nx',6),('ny',7),('nxy',8),('m1',9),('m2',10),('alpha_M',11),\
               ('mx_M',12),('my_M',13),('nx_M',14),('ny_M',15),  \
               ('sig_n',16),('sig_m',17),('eta_n',18),('eta_m',19), ( 'eta_tot', 16 ) ]
                
    font        = 'Courier 10'
    alignment   = 'right'
    format      = '%g'
   # column      = 'Elem.No.'
   # format      = '%d'
    index_text  = Property
#    index_image = Property
    
    def _get_index_text ( self ):
        return str( self.row )
        
    def _get_index_image ( self ):
        x, y, z = self.item
        if sqrt( (x - 0.5) ** 2 + (y - 0.5) ** 2 + (z - 0.5) ** 2 ) <= 0.25:
            return '@icons:red_ball'
            
        return None
    

class InfoShellFile( HasTraits ):
    '''
    Represent a single test specifying the design parameters.
    and access to the measured data.
    '''

    # thickness at border [m]
    D_b = Float( 0.040 )
    
    # thickness at top [m]
    D_t = Float( 0.030 )
    
    # tensile trength [MPa]
    f_ctk = Float( 1.6 )
    
    # flexural tensile trength [MPa]
    f_m = Float( 10.5 )
    
    data_file = Str

    result_arr_size = Int(4)
        
    NX = Array
    NY = Array
    MX = Array
    MY = Array

    @on_trait_change('data_file,result_arr_size')
    def _run_evaluation(self):
        
        # thickness at border [m]
        D_b = self.D_b
        
        # thickness at top [m]
        D_t = self.D_t
        
        # tensile trength [MPa]
        f_ctk = self.f_ctk
        
        # flexural tensile trength [MPa]
        f_m = self.f_m        
        # ------------------------------------------------------------
        # Read input txt-file and get data into numpy format
        # ------------------------------------------------------------
        # open input file containing raw data:
        file = open( self.data_file,'r')
        
        # open subsidary file for the conversion between "." and "," defined floats:
        a_data = open('adapted_data','w')
        
        # search and replace "," with ".":
        adapted = file.read().replace(",",".")
        
        # print 'adapted':
        a_data.write(adapted)
        a_data.close()
        
        # get subsidary file to retriev input data as array of floats:
        input_arr = loadtxt( 'adapted_data' , usecols=(-8,-7,-6,-5,-4,-3), delimiter='\t' )
        
        file = open( 'adapted_data','r')
        
        # moments
        mx  = input_arr[:,0]
        my  = input_arr[:,1]
        mxy = input_arr[:,2]
        
        
        # normalforces
        nx  = input_arr[:,3]
        ny  = input_arr[:,4]
        nxy = input_arr[:,5]
        
        
        # ------------------------------------------------------------
        # Get element number and node number
        # ------------------------------------------------------------
        elem_no = ones_like(mx)
        node_no  = ones_like(mx)
        
        cr_elem = 0
        for i in range( elem_no.shape[0] ):
            line = file.readline()
            if line.isspace():
                line = file.readline()
            value = split( line, '\t' )
            if value[0]:
                cr_elem = int(value[0])
                cr_node = int(value[1] )
            elem_no[i] = cr_elem
            node_no[i] = cr_node
            #if not line: break
        print "elem_no ",elem_no.shape
            
        row_no = arange(elem_no.shape[0])
        
        # ------------------------------------------------------------
        # Get area A, moment of inertia W
        # ------------------------------------------------------------
        # derive the element thickness based on the element number
        # (the division by 100 (=int) yields an integer)
        D_elem = D_b - (elem_no/100-1)/11 * (D_b - D_t)
        
        A = D_elem * 1.
        W = 1. * D_elem**2/6.
        
        
        # ------------------------------------------------------------
        # Index M: principle moments with corresponding normalforces
        # ------------------------------------------------------------
        
        # principle moments
        m1 = 0.5*(mx+my) + 0.5*sqrt((mx-my)**2 + 4*mxy**2)
        m2 = 0.5*(mx+my) - 0.5*sqrt((mx-my)**2 + 4*mxy**2)
        
        # principle angle of the moments:
        alpha_M = pi/2*ones_like(m1)
        alpha_M[m1 != m2] = arctan(mxy/(m2-mx))
        
        # transform the values in the principle direction:
        mx_M = 0.5*(mx+my) - 0.5*(mx-my)*cos(2*alpha_M) - mxy*sin(2*alpha_M)
        my_M = 0.5*(mx+my) + 0.5*(mx-my)*cos(2*alpha_M) + mxy*sin(2*alpha_M)
        
        nx_M = 0.5*(nx+ny) - 0.5*(nx-ny)*cos(2*alpha_M) - nxy*sin(2*alpha_M)
        ny_M = 0.5*(nx+ny) + 0.5*(nx-ny)*cos(2*alpha_M) + nxy*sin(2*alpha_M)
        
        
        # ------------------------------------------------------------
        # Index N: principle normal forces with corresponding moments
        # ------------------------------------------------------------
        
        # principle normal forces
        n1 = 0.5*(nx+ny) + 0.5*sqrt((nx-ny)**2 + 4*nxy**2)
        n2 = 0.5*(nx+ny) - 0.5*sqrt((nx-ny)**2 + 4*nxy**2)
        
        # principle angle of the moments:
        alpha_N = pi/2*ones_like(n1)
        alpha_N[n1 != n2] = arctan(nxy/(n2-nx))
        
        nx_N = 0.5*(nx+ny) - 0.5*(nx-ny)*cos(2*alpha_N) - nxy*sin(2*alpha_N)
        ny_N = 0.5*(nx+ny) + 0.5*(nx-ny)*cos(2*alpha_N) + nxy*sin(2*alpha_N)
        
        mx_N = 0.5*(mx+my) - 0.5*(mx-my)*cos(2*alpha_N) - mxy*sin(2*alpha_N)
        my_N = 0.5*(mx+my) + 0.5*(mx-my)*cos(2*alpha_N) + mxy*sin(2*alpha_N)
        
        
        # ------------------------------------------------------------
        # Evaluation of stresses
        # ------------------------------------------------------------ 
        # Index X: evaluation of stresses in transformed X-direction
        # Index Y: evaluation of stresses in transformed Y-direction
        
        # case: MX
        sig_n_MX = nx_M / A
        sig_m_MX = mx_M / W
        eta_n_MX = sig_n_MX / f_ctk
        eta_m_MX = sig_m_MX / f_m
        eta_tot_MX = eta_n_MX + eta_m_MX
        MX_ix = argsort( eta_tot_MX )
        
        # case: MY
        sig_n_MY = nx_M / A
        sig_m_MY = mx_M / W
        eta_n_MY = sig_n_MY / f_ctk
        eta_m_MY = sig_m_MY / f_m
        eta_tot_MY = eta_n_MY + eta_m_MY
        MY_ix = argsort( eta_tot_MY )
        
        # case: NX
        sig_n_NX = nx_N / A
        sig_m_NX = mx_N / W
        eta_n_NX = sig_n_NX / f_ctk
        eta_m_NX = sig_m_NX / f_m
        eta_tot_NX = eta_n_NX + eta_m_NX
        NX_ix = argsort( eta_tot_NX )
        
        # case: NY
        sig_n_NY = nx_N / A
        sig_m_NY = mx_N / W
        eta_n_NY = sig_n_NY / f_ctk
        eta_m_NY = sig_m_NY / f_m
        eta_tot_NY = eta_n_NY + eta_m_NY
        NY_ix = argsort( eta_tot_NY )
        
        
        NRES = self.result_arr_size
        self.NX = zeros( ( NRES, 21 ), dtype = 'float_' )
        self.NY = zeros( ( NRES, 21 ), dtype = 'float_' )
        self.MX = zeros( ( NRES, 21 ), dtype = 'float_' )
        self.MY = zeros( ( NRES, 21 ), dtype = 'float_' )

        # NX:
        self.NX[:,0] = row_no[ NX_ix[:NRES]] 
        self.NX[:,1] = elem_no[ NX_ix[:NRES]]
        self.NX[:,2] = node_no[ NX_ix[:NRES]]
        self.NX[:,3] = mx[ NX_ix[:NRES]]
        self.NX[:,4] = my[ NX_ix[:NRES]]
        self.NX[:,5] = mxy[ NX_ix[:NRES]]
        self.NX[:,6] = nx[ NX_ix[:NRES]]
        self.NX[:,7] = ny[ NX_ix[:NRES]]
        self.NX[:,8] = nxy[ NX_ix[:NRES]]
        #        
        self.NX[:,9] = n1[ NX_ix[:NRES]]
        self.NX[:,10] = n2[ NX_ix[:NRES]]
        self.NX[:,11] = alpha_N[ NX_ix[:NRES]]
        #
        self.NX[:,12] = mx_N[ NX_ix[:NRES]]
        self.NX[:,13] = my_N[ NX_ix[:NRES]]
        self.NX[:,14] = nx_N[ NX_ix[:NRES]]
        self.NX[:,15] = ny_N[ NX_ix[:NRES]]
        #
        self.NX[:,16] = sig_n_NX[ NX_ix[:NRES]]
        self.NX[:,17] = sig_m_NX[ NX_ix[:NRES]]
        self.NX[:,18] = eta_n_NX[ NX_ix[:NRES]]
        self.NX[:,19] = eta_m_NX[ NX_ix[:NRES]]
        self.NX[:,20] = eta_tot_NX[ NX_ix[:NRES]]
        
        # NY:
        self.NY[:,0] = row_no[ NY_ix[:NRES]] 
        self.NY[:,1] = elem_no[ NY_ix[:NRES]]
        self.NY[:,2] = node_no[ NY_ix[:NRES]]
        self.NY[:,3] = mx[ NY_ix[:NRES]]
        self.NY[:,4] = my[ NY_ix[:NRES]]
        self.NY[:,5] = mxy[ NY_ix[:NRES]]
        self.NY[:,6] = nx[ NY_ix[:NRES]]
        self.NY[:,7] = ny[ NY_ix[:NRES]]
        self.NY[:,8] = nxy[ NY_ix[:NRES]]
        #        
        self.NY[:,9] = n1[ NY_ix[:NRES]]
        self.NY[:,10] = n2[ NY_ix[:NRES]]
        self.NY[:,11] = alpha_N[ NY_ix[:NRES]]
        #
        self.NY[:,12] = mx_N[ NY_ix[:NRES]]
        self.NY[:,13] = my_N[ NY_ix[:NRES]]
        self.NY[:,14] = nx_N[ NY_ix[:NRES]]
        self.NY[:,15] = ny_N[ NY_ix[:NRES]]
        #
        self.NY[:,16] = sig_n_NY[ NY_ix[:NRES]]
        self.NY[:,17] = sig_m_NY[ NY_ix[:NRES]]
        self.NY[:,18] = eta_n_NY[ NY_ix[:NRES]]
        self.NY[:,19] = eta_m_NY[ NY_ix[:NRES]]
        self.NY[:,20] = eta_tot_NY[ NY_ix[:NRES]]
        
        # MX:
        self.MX[:,0] = row_no[ MX_ix[:NRES]] 
        self.MX[:,1] = elem_no[ MX_ix[:NRES]]
        self.MX[:,2] = node_no[ MX_ix[:NRES]]
        self.MX[:,3] = mx[ MX_ix[:NRES]]
        self.MX[:,4] = my[ MX_ix[:NRES]]
        self.MX[:,5] = mxy[ MX_ix[:NRES]]
        self.MX[:,6] = nx[ MX_ix[:NRES]]
        self.MX[:,7] = ny[ MX_ix[:NRES]]
        self.MX[:,8] = nxy[ MX_ix[:NRES]]
        #        
        self.MX[:,9] = m1[ MX_ix[:NRES]]
        self.MX[:,10] = m2[ MX_ix[:NRES]]
        self.MX[:,11] = alpha_M[ MX_ix[:NRES]]
        #
        self.MX[:,12] = mx_M[ MX_ix[:NRES]]
        self.MX[:,13] = my_M[ MX_ix[:NRES]]
        self.MX[:,14] = nx_M[ MX_ix[:NRES]]
        self.MX[:,15] = ny_M[ MX_ix[:NRES]]
        #
        self.MX[:,16] = sig_n_MX[ MX_ix[:NRES]]
        self.MX[:,17] = sig_m_MX[ MX_ix[:NRES]]
        self.MX[:,18] = eta_n_MX[ MX_ix[:NRES]]
        self.MX[:,19] = eta_m_MX[ MX_ix[:NRES]]
        self.MX[:,20] = eta_tot_MX[ MX_ix[:NRES]]
        
        # MY:
        self.MY[:,0] = row_no[ MY_ix[:NRES]] 
        self.MY[:,1] = elem_no[ MY_ix[:NRES]]
        self.MY[:,2] = node_no[ MY_ix[:NRES]]
        self.MY[:,3] = mx[ MY_ix[:NRES]]
        self.MY[:,4] = my[ MY_ix[:NRES]]
        self.MY[:,5] = mxy[ MY_ix[:NRES]]
        self.MY[:,6] = nx[ MY_ix[:NRES]]
        self.MY[:,7] = ny[ MY_ix[:NRES]]
        self.MY[:,8] = nxy[ MY_ix[:NRES]]
        #        
        self.MY[:,9] = m1[ MY_ix[:NRES]]
        self.MY[:,10] = m2[ MY_ix[:NRES]]
        self.MY[:,11] = alpha_M[ MY_ix[:NRES]]
        #
        self.MY[:,12] = mx_M[ MY_ix[:NRES]]
        self.MY[:,13] = my_M[ MY_ix[:NRES]]
        self.MY[:,14] = nx_M[ MY_ix[:NRES]]
        self.MY[:,15] = ny_M[ MY_ix[:NRES]]
        #
        self.MY[:,16] = sig_n_MY[ MY_ix[:NRES]]
        self.MY[:,17] = sig_m_MY[ MY_ix[:NRES]]
        self.MY[:,18] = eta_n_MY[ MY_ix[:NRES]]
        self.MY[:,19] = eta_m_MY[ MY_ix[:NRES]]
        self.MY[:,20] = eta_tot_MY[ MY_ix[:NRES]]        

                    
    traits_view = View(
                        Item( 'data_file', style = 'readonly' ),
                        Item( name = 'result_arr_size', label='row size'),
                        Tabbed(
                            Item( 'NX' , show_label = False, editor = TabularEditor( adapter = ArrayAdapterN() ) ),
                            Item( 'NY' , show_label = False, editor = TabularEditor( adapter = ArrayAdapterN() ) ),
                            Item( 'MX' , show_label = False, editor = TabularEditor( adapter = ArrayAdapterM() ) ),
                            Item( 'MY' , show_label = False, editor = TabularEditor( adapter = ArrayAdapterM() ) ),
                            scrollable = True,
                         ),
                      resizable = True,
                      scrollable = True,
                      height = 500,
                      width = 1000
                      )
   
isf = InfoShellFile( data_file = 'Rohdaten_GdG' )
#isf._run_evaluation()
isf.configure_traits()