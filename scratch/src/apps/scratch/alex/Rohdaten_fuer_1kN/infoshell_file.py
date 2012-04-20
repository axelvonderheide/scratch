
from enthought.traits.api import \
    HasTraits, Directory, List, Int, Float, Any, \
    on_trait_change, File, Constant, Instance, Trait, \
    Array, Str, Property, cached_property, WeakRef, \
    Dict, Button, Enum, Color
    
from enthought.util.home_directory import \
    get_home_directory

from enthought.traits.ui.api import \
    View, Item, DirectoryEditor, TabularEditor, HSplit, Tabbed, VGroup,\
    TableEditor, Group, ListEditor, VSplit, HSplit

from enthought.traits.ui.table_column import \
    ObjectColumn
    
from enthought.traits.ui.menu import \
    OKButton, CancelButton
    
from enthought.traits.ui.tabular_adapter \
    import TabularAdapter
    
from numpy import array, loadtxt, arange, sqrt, zeros, arctan, sin, cos, ones_like,\
                vstack, savetxt, hstack, argsort, fromstring, zeros_like, fix, copy  

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

class ArrayAdapter ( TabularAdapter ):
    
    columns = List()
    font        = 'Courier 10'
    alignment   = 'right'
    format      = '%5.2f'#'%g'
    even_bg_color = Color(0xE0E0FF)
    width = Float(80)
    
    # @temporary hack; is there a better solution?
    def get_format ( self, object, trait, row, column ):
        """ Returns the Python format string to use for a specified column.
        """
        if column == 0 or column == 1 or column == 2:
            a = self.format
            b = '%3d'
            self.format = b
            c = self._result_for( 'get_format', object, trait, row, column )
            self.format = a
            return c       
        else:
            return self._result_for( 'get_format', object, trait, row, column )    
    
class ArrayAdapterGM ( ArrayAdapter ):
    
    columns = [ ( 'Row.No.', 0), ( 'Elem.No.', 1 ), ( 'Node.No.', 2 ),
                ('mx', 3), ('my',4), ('mxy',5), ('nx',6), ('ny',7), ('nxy',8),\
                ('m1',9), ('m2',10), ('alpha_M',11),\
                ('mx_M',12), ('my_M',13), ('nx_M',14), ('ny_M',15),\
                ('sig_n',16), ('sig_m',17), ('eta_n',18), ('eta_m',19), ( 'eta_tot', 20 ) ]
                
class ArrayAdapterGN ( ArrayAdapter ):
    
    columns = [ ( 'Row.No.', 0 ), ( 'Elem.No.', 1 ), ( 'Node.No.', 2 ),\
                ('mx', 3), ('my',4), ('mxy',5), ('nx',6), ('ny',7), ('nxy',8),\
                ('n1',9), ('n2',10), ('alpha_N',11),\
                ('mx_N',12), ('my_N',13), ('nx_N',14), ('ny_N',15),  \
                ('sig_n',16), ('sig_m',17), ('eta_n',18), ('eta_m',19), ( 'eta_tot', 20 ) ]
    
class ArrayAdapterTM ( ArrayAdapter ):
    
    columns = [ ( 'Row.No.', 0), ( 'Elem.No.', 1 ), ( 'Node.No.', 2 ),\
                ('mx', 3), ('my',4), ('mxy',5), ('nx',6), ('ny',7), ('nxy',8),\
                ('m1',9), ('m2',10), ('alpha_M',11),\
                ('mx_M',12), ('my_M',13), ('nx_M',14), ('ny_M',15),  \
                ('e',16), ('m_Eds',17), ('f_t',18), ( 'f_Rtex', 19 ), ( 'n_tex', 20 ) ]
    
class ArrayAdapterTN ( ArrayAdapter ):
    
    columns = [ ( 'Row.No.', 0), ( 'Elem.No.', 1 ), ( 'Node.No.', 2 ), ('mx', 3), \
                ('my',4), ('mxy',5), ('nx',6), ('ny',7), ('nxy',8),\
                ('n1',9), ('n2',10), ('alpha_N',11),\
                ('mx_N',12), ('my_N',13), ('nx_N',14), ('ny_N',15),  \
                ('e',16), ('m_Eds',17), ('f_t',18), ( 'f_Rtex', 19 ), ( 'n_tex', 20 ) ]
       

class InfoShellFile( HasTraits ):
    '''
    Represent a single test specifying the design parameters.
    and access to the measured data.
    
    Made Assumptions for defined formulas:
    - reinforcement layers are spaced equally throughout the cross-section (d = 0,75*D_elem)
    - same reinforcement layers on top as on bottom of cross section (e = abs(m/n)
    '''

    # ------------------------------------------------------------
    # Geometry parameters:
    # ------------------------------------------------------------
    
    # thickness at border [m]
    D_b = Float( 0.040 )
    
    # thickness at top [m]
    D_t = Float( 0.030 )

    # ------------------------------------------------------------
    # GdG: Parameters for the uncracked state (G):
    # ------------------------------------------------------------
    
    # tensile strength [MPa]
    f_ctk = Float( 1.6 )

    # flexural tensile strength [MPa]
    f_m = Float( 10.5 )
    
    #--------------------------------------------------------
    # GdT: Prameters for the cracked state / reinforcement 
    #--------------------------------------------------------
    
    # INDEX l: longitudinal direction of the textile (MAG-02-02-06a)
    # tensile strength of the tensile specimen [N/mm2]
    # containing a gamma-factor of 1.5 and d long term reduction factor of 0.7
    f_td_l = 251
    
    # cross sectional area of the reinforcement [mm2/m]
    a_t_l = 71.65
    
    # INDEX q: orthogonal direction of the textile (MAG-02-02-06a)
    # tensile strength of the tensile specimen [kN/m]
    f_td_q = 238
    
    # cross sectional area of the reinforcement [mm2/m]
    a_t_q = 53.31
    
    # tensile strength of the textile reinforcement [kN/m]
    F_Rtex_l = Float( a_t_l * f_td_l / 1000.)  
    
    # tensile strength of the textile reinforcement [kN/m]
    F_Rtex_q = Float( a_t_q * f_td_q / 1000.)


    #------------------------------------------
    # Calculation of stresses and dimensioning:
    #------------------------------------------

    # raw input file
    data_file = Str

    # number of rows to be displayed
    result_arr_size = Int(4)
    
    # output arrays    
    NX_G = Array
    NY_G = Array
    MX_G = Array
    MY_G = Array

    NX_T = Array
    NY_T = Array
    MX_T = Array
    MY_T = Array

    ordering   = Enum("max_values", "row_no")
    
    eta_tot_max = Float()
    eta_tot_max_case = Str()
    n_tex_max = Float()
    n_tex_max_case = Str()
    
    @on_trait_change('data_file, result_arr_size, D_b, D_t, f_ctk, f_m, ordering')
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
        
        # search and replace "," with ".":
        adapted = file.read().replace(",",".")

        # open subsidary file for the conversion between "." and "," defined floats:
        a_data = open('adapted_data.txt','w')
        
        # print 'adapted':
        a_data.write(adapted)
        a_data.close()
        
        # get subsidary file to retriev input data as array of floats:
        input_arr = loadtxt( 'adapted_data.txt' , usecols=(-8,-7,-6,-5,-4,-3), delimiter='\t' )

#        r = np.loadtxt('record.dat', dtype={'names':('gender','age','weight'),
#                    'formats': ('S1','i4', 'f4')})


        # define arrays containing the information from the raw input file
        # moments:
        mx  = input_arr[:,0]
        my  = input_arr[:,1]
        mxy = input_arr[:,2]        
        # normalforces:
        nx  = input_arr[:,3]
        ny  = input_arr[:,4]
        nxy = input_arr[:,5]
                
        # ------------------------------------------------------------
        # Get element number and node number
        # ------------------------------------------------------------

        # open subsidary file to read element number and node numbers:
        file = open( 'adapted_data.txt','r')

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
        
        file.close()

        # remove the subsidary file from the current working directory
        cwd = os.getcwd()
        os_sep = os.sep
        os.remove(cwd + os_sep +'adapted_data.txt')
        
        
#        row_no = arange(elem_no.shape[0])
        
        # enumerate the rows as in infoCAD (also the empty rows have a number)
        # one element block in the adapted data file has 16 rows
        # first step derive the number of elements by devision by 16, then
        # create the same number of entries in the array 'aa' with step 17
        # starting with '1' 
        aa = arange(1,(elem_no.shape[0]/16)*17,17)
        bb = arange(16)
        cc = aa[:,None] + bb
        dd = cc.flatten()
        row_no = dd[:]
        
        # ------------------------------------------------------------
        # Index M: principle moments with corresponding normal forces
        # ------------------------------------------------------------
        
        # principle moments
        m1 = 0.5*(mx+my) + 0.5*sqrt((mx-my)**2 + 4*mxy**2)
        m2 = 0.5*(mx+my) - 0.5*sqrt((mx-my)**2 + 4*mxy**2)
        
        # principle angle of the moments:
        alpha_M = pi/2.*ones_like(m1)
        bool_M = m2 != mx
        alpha_M[bool_M] = arctan(mxy[bool_M]/(m2[bool_M]-mx[bool_M]))
#        for i in range(alpha_M.shape[0]):
#            if m2[i] != mx[i]:
#                alpha_M[i] = arctan(mxy[i]/(m2[i]-mx[i]))
        
        # transform the values in the principle direction:
        mx_M = 0.5*(my+mx) - 0.5*(my-mx)*cos(2*alpha_M) - mxy*sin(2*alpha_M)
        my_M = 0.5*(my+mx) + 0.5*(my-mx)*cos(2*alpha_M) + mxy*sin(2*alpha_M)
        
        nx_M = 0.5*(ny+nx) - 0.5*(ny-nx)*cos(2*alpha_M) - nxy*sin(2*alpha_M)
        ny_M = 0.5*(ny+nx) + 0.5*(ny-nx)*cos(2*alpha_M) + nxy*sin(2*alpha_M)
        
        
        # ------------------------------------------------------------
        # Index N: principle normal forces with corresponding moments
        # ------------------------------------------------------------
        
        # principle normal forces
        n1 = 0.5*(nx+ny) + 0.5*sqrt((nx-ny)**2 + 4*nxy**2)
        n2 = 0.5*(nx+ny) - 0.5*sqrt((nx-ny)**2 + 4*nxy**2)
        
        # principle angle of the moments:
        alpha_N = pi/2.*ones_like(n1)
        bool_N = n2 != nx
        alpha_N[bool_N] = arctan(nxy[bool_N]/(n2[bool_N]-nx[bool_N]))
        
        nx_N = 0.5*(ny+nx) - 0.5*(ny-nx)*cos(2*alpha_N) - nxy*sin(2*alpha_N)
        ny_N = 0.5*(ny+nx) + 0.5*(ny-nx)*cos(2*alpha_N) + nxy*sin(2*alpha_N)
        
        mx_N = 0.5*(my+mx) - 0.5*(my-mx)*cos(2*alpha_N) - mxy*sin(2*alpha_N)
        my_N = 0.5*(my+mx) + 0.5*(my-mx)*cos(2*alpha_N) + mxy*sin(2*alpha_N)
        
        
        # ------------------------------------------------------------
        # Index G: Evaluation of stresses
        # ------------------------------------------------------------ 
        # Index X: evaluation of stresses in transformed X-direction
        # Index Y: evaluation of stresses in transformed Y-direction

        # ------------------------------------------------------------
        # Get area A, moment of inertia W, calculate stresses
        # ------------------------------------------------------------
        # derive the element thickness based on the element number
        # (fix rounds to the next smallest integer)
        D_elem = D_b - (fix(elem_no/100)-1)/11 * (D_b - D_t)
        A = D_elem * 1. 
        W = 1. * D_elem**2/6. 
        
        # case: NX_G
        sig_n_NX = nx_N / A / 1000.
        sig_m_NX = abs( mx_N / W ) / 1000.
        eta_n_NX = sig_n_NX / f_ctk
        eta_m_NX = sig_m_NX / f_m
        eta_tot_NX = eta_n_NX + eta_m_NX

        # case: NY
        sig_n_NY = ny_N / A / 1000.
        sig_m_NY = abs( my_N / W ) / 1000.
        eta_n_NY = sig_n_NY / f_ctk
        eta_m_NY = sig_m_NY / f_m
        eta_tot_NY = eta_n_NY + eta_m_NY

        # case: MX_G stresses in [N/mm2]
        sig_n_MX = nx_M / A / 1000.
        sig_m_MX = abs( mx_M / W ) / 1000.
        eta_n_MX = sig_n_MX / f_ctk
        eta_m_MX = sig_m_MX / f_m
        eta_tot_MX = eta_n_MX + eta_m_MX
        
        # case: MY_G
        sig_n_MY = ny_M / A / 1000.
        sig_m_MY = abs( my_M / W ) / 1000.
        eta_n_MY = sig_n_MY / f_ctk
        eta_m_MY = sig_m_MY / f_m
        eta_tot_MY = eta_n_MY + eta_m_MY
        
        NRES = self.result_arr_size
        NX_G_ix = argsort( eta_tot_NX )[::-1]
        NY_G_ix = argsort( eta_tot_NY )[::-1]
        MX_G_ix = argsort( eta_tot_MX )[::-1]
        MY_G_ix = argsort( eta_tot_MY )[::-1]
    
        # calculate the maximum vale of eta_tot and the corresponding case 
        # (both are displayed in the view as read-only parameters) 
        eta_tot_NX_max = eta_tot_NX[ NX_G_ix[0]]
        eta_tot_NY_max = eta_tot_NY[ NY_G_ix[0]]
        eta_tot_MX_max = eta_tot_MX[ MX_G_ix[0]]
        eta_tot_MY_max = eta_tot_MY[ MY_G_ix[0]]
        #
        self.eta_tot_max = max( eta_tot_NX_max, eta_tot_NY_max, eta_tot_MX_max, eta_tot_MY_max )
        #
        if self.eta_tot_max == eta_tot_NX_max:
            self.eta_tot_max_case =  'G-NX' 
        elif self.eta_tot_max == eta_tot_NY_max:
            self.eta_tot_max_case =  'G-NY' 
        elif self.eta_tot_max == eta_tot_MX_max:
            self.eta_tot_max_case =  'G-MX' 
        elif self.eta_tot_max == eta_tot_MY_max:
            self.eta_tot_max_case =  'G-MY' 
        #
    
        if self.ordering == "row_no":
            print 'display all values (unordered)'
            NX_G_ix = arange(row_no.shape[0])
            NY_G_ix = arange(row_no.shape[0])
            MX_G_ix = arange(row_no.shape[0])
            MY_G_ix = arange(row_no.shape[0])
            NRES = row_no.shape[0]
        
        # ------------------------------------------------------------
        # G: prepare output array 
        # ------------------------------------------------------------

        self.NX_G = zeros( ( NRES, 21 ), dtype = 'float_' )
        self.NY_G = zeros( ( NRES, 21 ), dtype = 'float_' )
        self.MX_G = zeros( ( NRES, 21 ), dtype = 'float_' )
        self.MY_G = zeros( ( NRES, 21 ), dtype = 'float_' )

        # NX_G:
        self.NX_G[:,0] = row_no [ NX_G_ix[:NRES]] 
        self.NX_G[:,1] = elem_no[ NX_G_ix[:NRES]]
        self.NX_G[:,2] = node_no[ NX_G_ix[:NRES]]
        self.NX_G[:,3] = mx [ NX_G_ix[:NRES]]
        self.NX_G[:,4] = my [ NX_G_ix[:NRES]]
        self.NX_G[:,5] = mxy[ NX_G_ix[:NRES]]
        self.NX_G[:,6] = nx [ NX_G_ix[:NRES]]
        self.NX_G[:,7] = ny [ NX_G_ix[:NRES]]
        self.NX_G[:,8] = nxy[ NX_G_ix[:NRES]]
        self.NX_G[:,9] = n1 [ NX_G_ix[:NRES]]
        self.NX_G[:,10] = n2[ NX_G_ix[:NRES]]
        self.NX_G[:,11] = alpha_N[ NX_G_ix[:NRES]]
        self.NX_G[:,12] = mx_N[ NX_G_ix[:NRES]]
        self.NX_G[:,13] = my_N[ NX_G_ix[:NRES]]
        self.NX_G[:,14] = nx_N[ NX_G_ix[:NRES]]
        self.NX_G[:,15] = ny_N[ NX_G_ix[:NRES]]
        self.NX_G[:,16] = sig_n_NX[ NX_G_ix[:NRES]]
        self.NX_G[:,17] = sig_m_NX[ NX_G_ix[:NRES]]
        self.NX_G[:,18] = eta_n_NX[ NX_G_ix[:NRES]]
        self.NX_G[:,19] = eta_m_NX[ NX_G_ix[:NRES]]
        self.NX_G[:,20] = eta_tot_NX[ NX_G_ix[:NRES]]
        
        # NY_G:
        self.NY_G[:,0] = row_no [ NY_G_ix[:NRES]] 
        self.NY_G[:,1] = elem_no[ NY_G_ix[:NRES]]
        self.NY_G[:,2] = node_no[ NY_G_ix[:NRES]]
        self.NY_G[:,3] = mx [ NY_G_ix[:NRES]]
        self.NY_G[:,4] = my [ NY_G_ix[:NRES]]
        self.NY_G[:,5] = mxy[ NY_G_ix[:NRES]]
        self.NY_G[:,6] = nx [ NY_G_ix[:NRES]]
        self.NY_G[:,7] = ny [ NY_G_ix[:NRES]]
        self.NY_G[:,8] = nxy[ NY_G_ix[:NRES]]
        self.NY_G[:,9] = n1 [ NY_G_ix[:NRES]]
        self.NY_G[:,10] = n2[ NY_G_ix[:NRES]]
        self.NY_G[:,11] = alpha_N[ NY_G_ix[:NRES]]
        self.NY_G[:,12] = mx_N[ NY_G_ix[:NRES]]
        self.NY_G[:,13] = my_N[ NY_G_ix[:NRES]]
        self.NY_G[:,14] = nx_N[ NY_G_ix[:NRES]]
        self.NY_G[:,15] = ny_N[ NY_G_ix[:NRES]]
        self.NY_G[:,16] = sig_n_NY[ NY_G_ix[:NRES]]
        self.NY_G[:,17] = sig_m_NY[ NY_G_ix[:NRES]]
        self.NY_G[:,18] = eta_n_NY[ NY_G_ix[:NRES]]
        self.NY_G[:,19] = eta_m_NY[ NY_G_ix[:NRES]]
        self.NY_G[:,20] = eta_tot_NY[ NY_G_ix[:NRES]]
        
        # MX_G:
        self.MX_G[:,0] = row_no [ MX_G_ix[:NRES]] 
        self.MX_G[:,1] = elem_no[ MX_G_ix[:NRES]]
        self.MX_G[:,2] = node_no[ MX_G_ix[:NRES]]
        self.MX_G[:,3] = mx [ MX_G_ix[:NRES]]
        self.MX_G[:,4] = my [ MX_G_ix[:NRES]]
        self.MX_G[:,5] = mxy[ MX_G_ix[:NRES]]
        self.MX_G[:,6] = nx [ MX_G_ix[:NRES]]
        self.MX_G[:,7] = ny [ MX_G_ix[:NRES]]
        self.MX_G[:,8] = nxy[ MX_G_ix[:NRES]]
        self.MX_G[:,9] = m1 [ MX_G_ix[:NRES]]
        self.MX_G[:,10] = m2[ MX_G_ix[:NRES]]
        self.MX_G[:,11] = alpha_M[ MX_G_ix[:NRES]]
        self.MX_G[:,12] = mx_M[ MX_G_ix[:NRES]]
        self.MX_G[:,13] = my_M[ MX_G_ix[:NRES]]
        self.MX_G[:,14] = nx_M[ MX_G_ix[:NRES]]
        self.MX_G[:,15] = ny_M[ MX_G_ix[:NRES]]
        self.MX_G[:,16] = sig_n_MX[ MX_G_ix[:NRES]]
        self.MX_G[:,17] = sig_m_MX[ MX_G_ix[:NRES]]
        self.MX_G[:,18] = eta_n_MX[ MX_G_ix[:NRES]]
        self.MX_G[:,19] = eta_m_MX[ MX_G_ix[:NRES]]
        self.MX_G[:,20] = eta_tot_MX[ MX_G_ix[:NRES]]
        
        # MY_G:
        self.MY_G[:,0] = row_no [ MY_G_ix[:NRES]] 
        self.MY_G[:,1] = elem_no[ MY_G_ix[:NRES]]
        self.MY_G[:,2] = node_no[ MY_G_ix[:NRES]]
        self.MY_G[:,3] = mx [ MY_G_ix[:NRES]]
        self.MY_G[:,4] = my [ MY_G_ix[:NRES]]
        self.MY_G[:,5] = mxy[ MY_G_ix[:NRES]]
        self.MY_G[:,6] = nx [ MY_G_ix[:NRES]]
        self.MY_G[:,7] = ny [ MY_G_ix[:NRES]]
        self.MY_G[:,8] = nxy[ MY_G_ix[:NRES]]
        self.MY_G[:,9] = m1 [ MY_G_ix[:NRES]]
        self.MY_G[:,10] = m2[ MY_G_ix[:NRES]]
        self.MY_G[:,11] = alpha_M[ MY_G_ix[:NRES]]
        self.MY_G[:,12] = mx_M[ MY_G_ix[:NRES]]
        self.MY_G[:,13] = my_M[ MY_G_ix[:NRES]]
        self.MY_G[:,14] = nx_M[ MY_G_ix[:NRES]]
        self.MY_G[:,15] = ny_M[ MY_G_ix[:NRES]]
        self.MY_G[:,16] = sig_n_MY[ MY_G_ix[:NRES]]
        self.MY_G[:,17] = sig_m_MY[ MY_G_ix[:NRES]]
        self.MY_G[:,18] = eta_n_MY[ MY_G_ix[:NRES]]
        self.MY_G[:,19] = eta_m_MY[ MY_G_ix[:NRES]]
        self.MY_G[:,20] = eta_tot_MY[ MY_G_ix[:NRES]]        


        # ------------------------------------------------------------
        # Index T: dimensioning / number of reinforcements layers
        # ------------------------------------------------------------
        # Index X: dimensioning in transformed X-direction
        # Index Y: dimensioning in transformed Y-direction

        # ------------------------------------------------------------
        # Parameters for the cracked state (GdT):
        # ------------------------------------------------------------
        # element thickness is calculated above:
        # D_elem = D_b - (fix(elem_no/100)-1)/11 * (D_b - D_t)
        
        # (resultierende statische Nutzhoehe) 
        d = 0.75 * D_elem
        
        # (Abstand Schwereachse zur resultierende Bewehrungslage) 
        # chose the same amount of reinforcement at the top as at the bottom 
        # i.e. zs = zs1 = zs2
        zs = d - D_elem/2.
        
        # (Innerer Hebelarm) 
        z = 0.9 * d

        # ------------------------------------------------------------
        # case T: NX
        # ------------------------------------------------------------
        # (Exzentrizitaet)
        e_NX = abs( mx_N / nx_N )
        e_NX[nx_N == 0] = 0 # if normal force is zero set e to zero

        # moment at the height of the resulting reinforcement layer:
        m_Eds_NX = abs(mx_N) - zs * nx_N
        
        # tensile force in the reinforcement for bending and compression
        f_t_NX = m_Eds_NX / z + nx_N
        
        # check if the two conditions are true:
        cond1_NX = nx_N > 0
        cond2_NX = e_NX < zs
        bool_arr_NX = cond1_NX * cond2_NX
        # in case of pure tension in the cross section:
        f_t_NX[bool_arr_NX] = nx_N[bool_arr_NX] * (zs[bool_arr_NX] + e_NX[bool_arr_NX])/(zs[bool_arr_NX] + zs[bool_arr_NX]) 
        
        # angel of deflection of the textile reinforcement for dimensioning in x-direction
        # distinguished between longtudinal (l) and transversal (q) direction
        beta_l_N = pi/2 - abs(alpha_N)
        beta_q_N = abs(alpha_N)
        
        # resulting strength of the bi-directional textile considering the 
        # deflection of the reinforcement in the loading direction:
        f_Rtex_NX = self.F_Rtex_l * cos( beta_l_N) * (1 - beta_l_N / (pi/2) ) +\
                    self.F_Rtex_q * cos( beta_q_N) * (1 - beta_q_N / (pi/2) ) 
        
        # necessary number of reinfocement layers
        n_tex_NX = f_t_NX / f_Rtex_NX
        

        # ------------------------------------------------------------
        # case T: NY
        # ------------------------------------------------------------
        # (Exzentrizitaet)
        e_NY = abs( my_N / ny_N )
        e_NY[ny_N == 0] = 0 # if normal force is zero set e to zero
        
        # moment at the height of the resulting reinforcement layer:
        m_Eds_NY = abs(my_N) - zs * ny_N
        
        # tensile force in the reinforcement for bending and compression
        f_t_NY = m_Eds_NY / z + ny_N
        
        # check if the two conditions are true:
        cond1_NY = ny_N > 0
        cond2_NY = e_NY < zs
        bool_arr_NY = cond1_NY * cond2_NY
        # in case of pure tension in the cross section:
        f_t_NY[bool_arr_NY] = ny_N[bool_arr_NY] * (zs[bool_arr_NY] + e_NY[bool_arr_NY])/(zs[bool_arr_NY] + zs[bool_arr_NY]) 
        
        # angel of deflection of the textile reinforcement for dimensioning in x-direction
        # distinguished between longtudinal (l) and transversal (q) direction
        beta_l_N = abs(alpha_N)
        beta_q_N = pi/2 - abs(alpha_N)
        
        # resulting strength of the bi-directional textile considering the 
        # deflection of the reinforcement in the loading direction:
        f_Rtex_NY = self.F_Rtex_l * cos( beta_l_N) * (1 - beta_l_N / (pi/2) ) +\
                    self.F_Rtex_q * cos( beta_q_N) * (1 - beta_q_N / (pi/2) )         

        # necessary number of reinfocement layers
        n_tex_NY = f_t_NY / f_Rtex_NY

        
        # ------------------------------------------------------------
        # case T: MX
        # ------------------------------------------------------------
        # (Exzentrizitaet)
        e_MX = abs( mx_M / nx_M )
        e_MX[nx_M == 0] = 0 # if normal force is zero set e to zero
        
        # moment at the height of the resulting reinforcement layer:
        m_Eds_MX = abs(mx_M) - zs * nx_M
        
        # tensile force in the reinforcement for bending and compression
        f_t_MX = m_Eds_MX / z + nx_M
        
        # check if the two conditions are true:
        cond1_MX = nx_M > 0
        cond2_MX = e_MX < zs
        bool_arr_MX = cond1_MX * cond2_MX
        # in case of pure tension in the cross section:
        f_t_MX[bool_arr_MX] = nx_M[bool_arr_MX] * (zs[bool_arr_MX] + e_MX[bool_arr_MX])/(zs[bool_arr_MX] + zs[bool_arr_MX]) 
        
        # angel of deflection of the textile reinforcement for dimensioning in x-direction
        # distinguished between longtudinal (l) and transversal (q) direction
        beta_l_N = pi/2 - abs(alpha_M)
        beta_q_N = abs(alpha_M)
        
        # resulting strength of the bi-directional textile considering the 
        # deflection of the reinforcement in the loading direction:
        f_Rtex_MX = self.F_Rtex_l * cos( beta_l_N) * (1 - beta_l_N / (pi/2) ) +\
                    self.F_Rtex_q * cos( beta_q_N) * (1 - beta_q_N / (pi/2) ) 
 
        # necessary number of reinfocement layers
        n_tex_MX = f_t_MX / f_Rtex_MX
        
        
        # ------------------------------------------------------------
        # case T: MY
        # ------------------------------------------------------------
        # (Exzentrizitaet)
        e_MY = abs( my_M / ny_M )
        e_MY[ny_M == 0] = 0 # if normal force is zero set e to zero
   
        # moment at the height of the resulting reinforcement layer:
        m_Eds_MY = abs(my_M) - zs * ny_M
        
        # tensile force in the reinforcement for bending and compression
        f_t_MY = m_Eds_MY / z + ny_M
        
        # check if the two conditions are true:
        cond1_MY = ny_M > 0
        cond2_MY = e_MY < zs
        bool_arr_MY = cond1_MY * cond2_MY
        # in case of pure tension in the cross section:
        f_t_MY[bool_arr_MY] = ny_M[bool_arr_MY] * (zs[bool_arr_MY] + e_MY[bool_arr_MY])/(zs[bool_arr_MY] + zs[bool_arr_MY]) 
        
        # angel of deflection of the textile reinforcement for dimensioning in x-direction
        # distinguished between longtudinal (l) and transversal (q) direction
        beta_l_N = abs(alpha_M)
        beta_q_N = pi/2 - abs(alpha_M)
        
        # resulting strength of the bi-directional textile considering the 
        # deflection of the reinforcement in the loading direction:
        f_Rtex_MY = self.F_Rtex_l * cos( beta_l_N) * (1 - beta_l_N / (pi/2) ) +\
                    self.F_Rtex_q * cos( beta_q_N) * (1 - beta_q_N / (pi/2) )         
        
        # necessary number of reinfocement layers
        n_tex_MY = f_t_MY / f_Rtex_MY
        
        
        # ------------------------------------------------------------
        # T: prepare output array 
        # ------------------------------------------------------------

        # sort the values depending on the maximum numer of reinfocement layers
        NX_T_ix = argsort( n_tex_NX )[::-1]
        NY_T_ix = argsort( n_tex_NY )[::-1]
        MX_T_ix = argsort( n_tex_MX )[::-1]
        MY_T_ix = argsort( n_tex_MY )[::-1]

        # calculate the maximum vale of eta_tot and the corresponding case 
        # (both are displayed in the view as read-only parameters) 
        n_tex_NX_max = n_tex_NX[ NX_T_ix[0]]
        n_tex_NY_max = n_tex_NY[ NY_T_ix[0]]
        n_tex_MX_max = n_tex_MX[ MX_T_ix[0]]
        n_tex_MY_max = n_tex_MY[ MY_T_ix[0]]
        #
        self.n_tex_max = max( n_tex_NX_max, n_tex_NY_max, n_tex_MX_max, n_tex_MY_max )
        #
        if self.n_tex_max == n_tex_NX_max:
            self.n_tex_max_case =  'T-NX' 
        elif self.n_tex_max == n_tex_NY_max:
            self.n_tex_max_case =  'T-NY' 
        elif self.n_tex_max == n_tex_MX_max:
            self.n_tex_max_case =  'T-MX' 
        elif self.n_tex_max == n_tex_MY_max:
            self.n_tex_max_case =  'T-MY' 
        #
        elif self.ordering == "row_no":
            NX_T_ix = arange(row_no.shape[0])
            NY_T_ix = arange(row_no.shape[0])
            MX_T_ix = arange(row_no.shape[0])
            MY_T_ix = arange(row_no.shape[0])
        #


        self.NX_T = zeros( ( NRES, 22 ), dtype = 'float_' )
        self.NY_T = zeros( ( NRES, 22 ), dtype = 'float_' )
        self.MX_T = zeros( ( NRES, 22 ), dtype = 'float_' )
        self.MY_T = zeros( ( NRES, 22 ), dtype = 'float_' )
        
        # NX_T:
        self.NX_T[:,0] = row_no [ NX_T_ix[:NRES]] 
        self.NX_T[:,1] = elem_no[ NX_T_ix[:NRES]]
        self.NX_T[:,2] = node_no[ NX_T_ix[:NRES]]
        self.NX_T[:,3] = mx [ NX_T_ix[:NRES]]
        self.NX_T[:,4] = my [ NX_T_ix[:NRES]]
        self.NX_T[:,5] = mxy[ NX_T_ix[:NRES]]
        self.NX_T[:,6] = nx [ NX_T_ix[:NRES]]
        self.NX_T[:,7] = ny [ NX_T_ix[:NRES]]
        self.NX_T[:,8] = nxy[ NX_T_ix[:NRES]]
        self.NX_T[:,9] = n1 [ NX_T_ix[:NRES]]
        self.NX_T[:,10] = n2[ NX_T_ix[:NRES]]
        self.NX_T[:,11] = alpha_N[ NX_T_ix[:NRES]]
        self.NX_T[:,12] = mx_N[ NX_T_ix[:NRES]]
        self.NX_T[:,13] = my_N[ NX_T_ix[:NRES]]
        self.NX_T[:,14] = nx_N[ NX_T_ix[:NRES]]
        self.NX_T[:,15] = ny_N[ NX_T_ix[:NRES]]
        self.NX_T[:,16] = e_NX[ NX_T_ix[:NRES]]
        self.NX_T[:,17] = m_Eds_NX[ NX_T_ix[:NRES]]
        self.NX_T[:,18] = f_t_NX  [ NX_T_ix[:NRES]]
        self.NX_T[:,19] = f_Rtex_NX[ NX_T_ix[:NRES]]
        self.NX_T[:,20] = n_tex_NX [ NX_T_ix[:NRES]]
        
        # NY_T:
        self.NY_T[:,0] = row_no [ NY_T_ix[:NRES]] 
        self.NY_T[:,1] = elem_no[ NY_T_ix[:NRES]]
        self.NY_T[:,2] = node_no[ NY_T_ix[:NRES]]
        self.NY_T[:,3] = mx [ NY_T_ix[:NRES]]
        self.NY_T[:,4] = my [ NY_T_ix[:NRES]]
        self.NY_T[:,5] = mxy[ NY_T_ix[:NRES]]
        self.NY_T[:,6] = nx [ NY_T_ix[:NRES]]
        self.NY_T[:,7] = ny [ NY_T_ix[:NRES]]
        self.NY_T[:,8] = nxy[ NY_T_ix[:NRES]]
        self.NY_T[:,9] = n1 [ NY_T_ix[:NRES]]
        self.NY_T[:,10] = n2[ NY_T_ix[:NRES]]
        self.NY_T[:,11] = alpha_N[ NY_T_ix[:NRES]]
        self.NY_T[:,12] = mx_N[ NY_T_ix[:NRES]]
        self.NY_T[:,13] = my_N[ NY_T_ix[:NRES]]
        self.NY_T[:,14] = nx_N[ NY_T_ix[:NRES]]
        self.NY_T[:,15] = ny_N[ NY_T_ix[:NRES]]
        self.NY_T[:,16] = e_NY[ NY_T_ix[:NRES]]
        self.NY_T[:,17] = m_Eds_NY[ NY_T_ix[:NRES]]
        self.NY_T[:,18] = f_t_NY  [ NY_T_ix[:NRES]]
        self.NY_T[:,19] = f_Rtex_NY[ NY_T_ix[:NRES]]
        self.NY_T[:,20] = n_tex_NY [ NY_T_ix[:NRES]]
        
        # MX_T:
        self.MX_T[:,0] = row_no [ MX_T_ix[:NRES]] 
        self.MX_T[:,1] = elem_no[ MX_T_ix[:NRES]]
        self.MX_T[:,2] = node_no[ MX_T_ix[:NRES]]
        self.MX_T[:,3] = mx [ MX_T_ix[:NRES]]
        self.MX_T[:,4] = my [ MX_T_ix[:NRES]]
        self.MX_T[:,5] = mxy[ MX_T_ix[:NRES]]
        self.MX_T[:,6] = nx [ MX_T_ix[:NRES]]
        self.MX_T[:,7] = ny [ MX_T_ix[:NRES]]
        self.MX_T[:,8] = nxy[ MX_T_ix[:NRES]]
        self.MX_T[:,9] = m1 [ MX_T_ix[:NRES]]
        self.MX_T[:,10] = m2[ MX_T_ix[:NRES]]
        self.MX_T[:,11] = alpha_M[ MX_T_ix[:NRES]]
        self.MX_T[:,12] = mx_M[ MX_T_ix[:NRES]]
        self.MX_T[:,13] = my_M[ MX_T_ix[:NRES]]
        self.MX_T[:,14] = nx_M[ MX_T_ix[:NRES]]
        self.MX_T[:,15] = ny_M[ MX_T_ix[:NRES]]
        self.MX_T[:,16] = e_MX[ MX_T_ix[:NRES]]
        self.MX_T[:,17] = m_Eds_MX[ MX_T_ix[:NRES]]
        self.MX_T[:,18] = f_t_MX  [ MX_T_ix[:NRES]]
        self.MX_T[:,19] = f_Rtex_MX[ MX_T_ix[:NRES]]
        self.MX_T[:,20] = n_tex_MX [ MX_T_ix[:NRES]]
       
        # MY_T:
        self.MY_T[:,0] = row_no [ MY_T_ix[:NRES]] 
        self.MY_T[:,1] = elem_no[ MY_T_ix[:NRES]]
        self.MY_T[:,2] = node_no[ MY_T_ix[:NRES]]
        self.MY_T[:,3] = mx [ MY_T_ix[:NRES]]
        self.MY_T[:,4] = my [ MY_T_ix[:NRES]]
        self.MY_T[:,5] = mxy[ MY_T_ix[:NRES]]
        self.MY_T[:,6] = nx [ MY_T_ix[:NRES]]
        self.MY_T[:,7] = ny [ MY_T_ix[:NRES]]
        self.MY_T[:,8] = nxy[ MY_T_ix[:NRES]]
        self.MY_T[:,9] = m1 [ MY_T_ix[:NRES]]
        self.MY_T[:,10] = m2[ MY_T_ix[:NRES]]
        self.MY_T[:,11] = alpha_M[ MY_T_ix[:NRES]]
        self.MY_T[:,12] = mx_M[ MY_T_ix[:NRES]]
        self.MY_T[:,13] = my_M[ MY_T_ix[:NRES]]
        self.MY_T[:,14] = nx_M[ MY_T_ix[:NRES]]
        self.MY_T[:,15] = ny_M[ MY_T_ix[:NRES]]
        self.MY_T[:,16] = e_MY[ MY_T_ix[:NRES]]
        self.MY_T[:,17] = m_Eds_MY [ MY_T_ix[:NRES]]
        self.MY_T[:,18] = f_t_MY   [ MY_T_ix[:NRES]]
        self.MY_T[:,19] = f_Rtex_MY[ MY_T_ix[:NRES]]
        self.MY_T[:,20] = n_tex_MY [ MY_T_ix[:NRES]]
    
    
    traits_view = View( Item( 'data_file', label='Evaluated input file ', style = 'readonly' ),
                        Item( name = 'D_b', label='Thickness at border [m]: D_b        '),
                        Item( name = 'D_t', label='Thickness at top [m]: D_t         '), 
                        Item( name = 'f_ctk', label='Tensile strength [MPa]:  f_ctk      '),
                        Item( name = 'f_m', label='Flexural tensile trength [MPa]:  f_m        '),
                        Item( name = 'F_Rtex_l', label='Strength textil (longitudinal) [kN/m]:  F_Rtex_l '),
                        Item( name = 'F_Rtex_q', label='Strength textil (transversal) [kN/m]: F_Rtex_q '),
                        HSplit(Item( name = 'eta_tot_max', label="Maximum value of 'eta_tot_max'                     ", style = 'readonly' ),
                               Item( name = 'eta_tot_max_case', label = 'reached in case', style = 'readonly' )),
                        HSplit(Item( name = 'n_tex_max', label="Maximum value of 'n_tex_max'                       ", style = 'readonly' ),
                               Item( name = 'n_tex_max_case', label = 'reached in case', style = 'readonly' ) ),
                        HSplit(Item('ordering', label='Use ordering', style = 'custom' ),
                               Item( name = 'result_arr_size', label='Numer of selected rows'),),
                        Tabbed(
                            Item( 'NX_G' , label = "G-NX", show_label = False, editor = TabularEditor( adapter = ArrayAdapterGN() ) ),
                            Item( 'NY_G' , label = "G-NY", show_label = False, editor = TabularEditor( adapter = ArrayAdapterGN() ) ),
                            Item( 'MX_G' , label = "G-MX", show_label = False, editor = TabularEditor( adapter = ArrayAdapterGM() ) ),
                            Item( 'MY_G' , label = "G-MY", show_label = False, editor = TabularEditor( adapter = ArrayAdapterGM() ) ),
                            Item( 'NX_T' , label = "T-NX", show_label = False, editor = TabularEditor( adapter = ArrayAdapterTN() ) ),
                            Item( 'NY_T' , label = "T-NY", show_label = False, editor = TabularEditor( adapter = ArrayAdapterTN() ) ),
                            Item( 'MX_T' , label = "T-MX", show_label = False, editor = TabularEditor( adapter = ArrayAdapterTM() ) ),
                            Item( 'MY_T' , label = "T-MY", show_label = False, editor = TabularEditor( adapter = ArrayAdapterTM() ) ),
                            scrollable = False,
                         ),
                      resizable = True,
                      scrollable = True,
                      height = 1000,
                      width = 1100
                      )
   
   
if __name__ == '__main__':
    isf = InfoShellFile( data_file = 'Ap_Asymmetrisch_Feld_links_1kN.rtf' )
    #isf._run_evaluation()
    isf.configure_traits()
