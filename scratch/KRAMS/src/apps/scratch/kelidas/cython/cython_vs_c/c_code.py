'''
Created on Nov 2, 2011

@author: kelidas
'''


from scipy.weave import inline, converters
import platform
import os


def spirrid_c(arg_list, c_params):

    C_code = '''
        for( int i_eps = 0; i_eps < n_eps; i_eps++){
            double eps = *( e_arr + i_eps );
            double mu_q(0);
            double q(0);
            for( int i_lambd = 0; i_lambd < n_int; i_lambd++){
                double lambd = *( Theta_la + i_lambd );
                for( int i_xi = 0; i_xi < n_int; i_xi++){
                    double xi = *( Theta_xi + i_xi );  
                        q = ( eps ) /( 1 + lambd );          
                        double eps_ = xi - q;
                        // Computation of the q( ... ) function
                        if ( eps_ < 0 || eps_ > xi ){
                        q = 0.0;
                        }
                // Store the values in the grid
                mu_q +=  q * *( g_la_pdf + i_lambd)* *( g_xi_pdf + i_xi);
            };
            };
        *(mu_q_arr + i_eps) = mu_q * d_la * d_xi;
        };
        '''

    compiler = 'gcc'

    compiler_verbose = 0

    conv = converters.default

    if platform.system() == 'Linux':
        #os.environ['CC'] = 'gcc-4.1'
        #os.environ['CXX'] = 'g++-4.1'
        os.environ['OPT'] = '-DNDEBUG -g -fwrapv -O3'
    elif platform.system() == 'Windows':
        # it is not Linux - just let it go and suffer
        pass

    return inline(C_code, arg_list, local_dict = c_params,
                                       type_converters = conv, compiler = compiler,
                                       verbose = compiler_verbose)
