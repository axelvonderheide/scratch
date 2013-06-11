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
# Created on May 25, 2010 by: rch

# @todo: the order of parameters must be defined - the only possibility is
# to define a list in the response function as it was the case orginally. 
# The current implementation # does not allow for the call to array-
# based RF evaluation.
#
# @todo: Compiled call to rf_grid for the calculation of the standard deviation.
#

from enthought.traits.api import \
    HasTraits, Instance, List, Property, Array, Int, Any, cached_property, Dict, \
    Event, on_trait_change, Bool, Float, WeakRef, Str

from numpy import \
    ogrid, trapz, linspace, array, sum, arange, zeros_like, \
    ones_like, max, argmax, sqrt, ones, copy as ncopy

import copy

from stats.spirrid.i_rf import \
    IRF

from stats.spirrid.rf_filament import \
    Filament

from scipy.weave import \
    inline, converters

from stats.pdistrib.pdistrib import \
    PDistrib

from string import \
    split

from types import \
    ListType

import os

import time

from mathkit.mfn.mfn_line.mfn_line import \
    MFnLineArray

def orthogonalize( arr_list ):
    '''Orthogonalize a list of one-dimensional arrays.
    '''
    n_arr = len( arr_list )
    ogrid = []
    for i, arr in enumerate( arr_list ):
        shape = ones( ( n_arr, ), dtype='int' )
        shape[i] = len( arr )
        arr_i = ncopy( arr ).reshape( tuple( shape ) )
        ogrid.append( arr_i )
    return ogrid

class RV( HasTraits ):
    '''Class representing the definition and discretization of a random variable.
    '''
    name = Str

    pd = Instance( PDistrib )

    n_int = Int( 30 )

    # index within the randomization
    idx = Int( 0 )

    theta_arr = Property( Array( 'float_' ), depends_on='distr' )
    @cached_property
    def _get_theta_arr( self ):
        return self.pd.x_array

    pdf_arr = Property( Array( 'float_' ), depends_on='distr' )
    @cached_property
    def _get_pdf_arr( self ):
        return self.pd.pdf_array

    pdf_theta_arr = Property( Array( 'float_' ), depends_on='distr' )
    @cached_property
    def _get_pdf_theta_arr( self ):
        pdf_theta_arr = self.pdf_arr * ( self.theta_arr[1] - self.theta_arr[0] )
        pdf_theta_arr[ 0] *= .5
        pdf_theta_arr[-1] *= .5
        return pdf_theta_arr

class SPIRRID( HasTraits ):
    '''Multidimensional statistical integration.
    
    Its name SPIRRID is an acronym for 
    Set of Parallel Independent Random Responses with Identical Distributions
    
    The package implements the evaluation of an integral over a set of 
    random variables affecting a response function RF and distributed 
    according to a probabilistic distribution PDistrib.
    
    The input parameters are devided in four categories in order
    to define state consistency of the evaluation. The outputs 
    are define as cached properties that are reevaluated in response
    to changes in the inputs.
    
    The following events accummulate changes in the input parameters of spirrid:
    rf_change - change in the response function
    rand_change - change in the randomization
    conf_change - change in the configuration of the algorithm
    eps_change - change in the studied range of the process control variable       
    '''
    #--------------------------------------------------------------------
    # Response function 
    #--------------------------------------------------------------------
    # @todo define the IRF - interface
    #
    rf = Instance( HasTraits )

    #--------------------------------------------------------------------
    # Specification of random parameters 
    #--------------------------------------------------------------------
    # 
    rv_dict = Dict
    def add_rv( self, variable, distribution='uniform',
                loc=0., scale=1., shape=1., n_int=30 ):
        '''Declare a variable as random 
        '''
        if variable not in self.rf.param_keys:
            raise AssertionError, 'parameter %s not defined by the response function' \
                % variable

        params_with_distr = self.rf.traits( distr=lambda x: type( x ) == ListType
                                            and distribution in x )
        if variable not in params_with_distr:
            raise AssertionError, 'distribution type %s not allowed for parameter %s' \
                % ( distribution, variable )

        # @todo - let the RV take care of PDistrib specification.
        # isolate the dirty two-step definition of the distrib from spirrid 
        #
        pd = PDistrib( distr_choice=distribution, n_segments=n_int )
        pd.distr_type.set( scale=scale, shape=shape, loc=loc )
        self.rv_dict[variable] = RV( name=variable, pd=pd, n_int=n_int )

    def del_rv( self, variable ):
        '''Delete declaration of random variable
        '''
        del self.rv_dict[ variable ]

    # subsidiary methods for sorted access to the random variables.
    # (note dictionary has not defined order of its items)
    rv_keys = Property( List, depends_on='rv_dict' )
    @cached_property
    def _get_rv_keys( self ):
        rv_keys = sorted( self.rv_dict.keys() )
        # the random variable gets an index based on the 
        # sorted keys
        for idx, key in enumerate( rv_keys ):
            self.rv_dict[ key ].idx = idx
        return rv_keys

    rv_list = Property( List, depends_on='rv_dict' )
    @cached_property
    def _get_rv_list( self ):
        return map( self.rv_dict.get, self.rv_keys )

    #---------------------------------------------------------------------------------------------
    # Range of the control process variable epsilon
    # Define particular control variable points with
    # the cv array or an equidistant range with (min, max, n)
    #---------------------------------------------------------------------------------------------
    
    cv = Array( eps_range=True )
    min_eps = Float( 0.0, eps_range=True )
    max_eps = Float( 0.0, eps_range=True )
    n_eps = Float( 80, eps_range=True )

    #--------------------------------------------------------------------
    # Define which changes in the response function and in the 
    # statistical parameters are relevant for reevaluation of the response
    #--------------------------------------------------------------------
    rf_change = Event
    @on_trait_change( 'rf.+distr' )
    def _set_rf_change( self ):
        self.rf_change = True

    rand_change = Event
    @on_trait_change( 'rand_params, rv_dict' )
    def _set_rand_change( self ):
        self.rand_change = True

    conf_change = Event
    @on_trait_change( '+alg_option' )
    def _set_conf_change( self ):
        self.conf_change = True

    eps_change = Event
    @on_trait_change( '+eps_range' )
    def _set_eps_change( self ):
        self.eps_change = True

    # Dictionary with key = rf parameters
    # and values = default param values for the resp func 
    #
    param_dict = Property( Dict, depends_on='rf_change, rand_change' )
    @cached_property
    def _get_param_dict( self ):
        '''Gather all the traits with the metadata distr specified.
        '''
        dict = {}
        for name, value in zip( self.rf.param_keys, self.rf.param_list ):
            rv = self.rv_dict.get( name, None )
            if rv == None:
                dict[ name ] = value
            else:
                dict[ name ] = self.theta_ogrid[ rv.idx ]
        return dict

    # Constant parameters
    #
    const_param_dict = Property( Dict, depends_on='rf_change, rand_change' )
    @cached_property
    def _get_const_param_dict( self ):
        const_param_dict = {}
        for name, v in zip( self.rf.param_keys, self.rf.param_list ):
            if name not in self.rv_keys:
                const_param_dict[ name ] = v
        return const_param_dict

    # Discretized statistical domain
    # 
    theta_ogrid = Property( depends_on='rf_change, rand_change' )
    @cached_property
    def _get_theta_ogrid( self ):
        '''Get orthogonal list of arrays with discretized RVs.
        '''
        theta_arr_list = [ rv.theta_arr for rv in self.rv_list ]
        return orthogonalize( theta_arr_list )

    #---------------------------------------------------------------------------------
    # PDF arrays oriented in enumerated dimensions - broadcasting possible
    #---------------------------------------------------------------------------------
    pdf_ogrid = Property( depends_on='rf_change, rand_change' )
    @cached_property
    def _get_pdf_ogrid( self ):
        '''Get orthogonal list of arrays with PDF values of RVs.
        '''
        pdf_arr_list = [ rv.pdf_arr for rv in self.rv_list ]
        return orthogonalize( pdf_arr_list )

    #---------------------------------------------------------------------------------
    # PDF * Theta arrays oriented in enumerated dimensions - broadcasting possible
    #---------------------------------------------------------------------------------
    pdf_theta_ogrid = Property( depends_on='rf_change, rand_change' )
    @cached_property
    def _get_pdf_theta_ogrid( self ):
        '''Get orthogonal list of arrays with PDF * Theta product of.
        '''
        pdf_theta_arr_list = [ rv.pdf_theta_arr for rv in self.rv_list ]
        return orthogonalize( pdf_theta_arr_list )
  
    #---------------------------------------------------------------------------------
    # PDF grid - mutually multiplied arrays of PDF
    #---------------------------------------------------------------------------------
    pdf_theta_grid = Property( depends_on='rf_change, rand_change' )
    @cached_property
    def _get_pdf_theta_grid( self ):
        return reduce( lambda x, y: x * y, self.pdf_theta_ogrid )
     
    #------------------------------------------------------------------------------------
    # Configuration of the algorithm
    #------------------------------------------------------------------------------------
    # 
    # cached_pdf_theta_grid:
    # If set to True, the cross product between the pdf values of all random variables
    # will be precalculated and stored in an n-dimensional grid
    # otherwise the product is performed for every epsilon in the inner loop anew
    # 
    cached_qg = Bool( True, alg_option=True )

    # compiled_eps_loop:
    # If set True, the loop over the control variable epsilon is compiled
    # otherwise, python loop is used.
    compiled_eps_loop = Bool( True, alg_option=True )

    # compiled_qg_loop:
    # If set True, the integration loop over the product between the response function
    # and the pdf . theta product is performed in c
    # otherwise the numpy arrays are used.
    compiled_qg_loop = Bool( True, alg_option=True )

    arg_list = Property( depends_on='rf_change, rand_change, conf_change' )
    @cached_property
    def _get_arg_list( self ):

        arg_list = []
        # create argument string for inline function
        if self.compiled_eps_loop:
            arg_list += [ 'mu_q_arr', 'e_arr' ]
        else:
            arg_list.append( 'e' )

        arg_list += ['%s_flat' % name for name in self.rv_keys ]

        if self.cached_qg:
            arg_list += [ 'pdf_theta_grid' ]
        else:
            arg_list += [ '%s_pdf' % name for name in self.rv_keys ]

        return arg_list

    C_code_qg = Property( depends_on='rf_change, rand_change, conf_change' )
    @cached_property
    def _get_C_code_qg( self ):
        if self.cached_qg: # q_g - blitz matrix used to store the grid
            code_str = '\tdouble pdf = pdf_theta_grid(' + \
                       ','.join( [ 'i_%s' % name
                                  for name in self.rv_keys ] ) + \
                       ');\n'
        else: # qg
            code_str = '\tdouble pdf = ' + \
                       '*'.join( [ ' *( %s_pdf + i_%s)' % ( name, name )
                                  for name in self.rv_keys ] ) + \
                       ';\n'
        return code_str

    #------------------------------------------------------------------------------------
    # Configurable generation of C-code for mean curve evaluation
    #------------------------------------------------------------------------------------
    C_code = Property( depends_on='rf_change, rand_change, conf_change, eps_change' )
    @cached_property
    def _get_C_code( self ):

        code_str = ''
        if self.compiled_eps_loop:

            # create code string for inline function
            #
            code_str += 'for( int i_eps = 0; i_eps < %g; i_eps++){\n' % self.n_eps

            if self.cached_qg:

                # multidimensional index needed for pdf_theta_grid 
                # use blitz arrays must be used also for other arrays
                #
                code_str += 'double eps = e_arr( i_eps );\n'

            else:
                # pointer access possible for single dimensional arrays 
                # use the pointer arithmetics for accessing the pdfs
                code_str += '\tdouble eps = *( e_arr + i_eps );\n'

        else:

            # create code string for inline function
            #
            code_str += 'double eps = e;\n'

        code_str += 'double mu_q(0);\n'
        code_str += 'double q(0);\n'

        code_str += '#line 100\n'
        # create code for constant params
        for name, value in self.const_param_dict.items():
            code_str += 'double %s = %g;\n' % ( name, value )

        # generate loops over random params

        for rv in self.rv_list:

            name = rv.name
            n_int = rv.n_int

            # create the loop over the random variable
            #
            code_str += 'for( int i_%s = 0; i_%s < %g; i_%s++){\n' % ( name, name, n_int, name )
            if self.cached_qg:

                # multidimensional index needed for pdf_grid - use blitz arrays
                #
                code_str += '\tdouble %s = %s_flat( i_%s );\n' % ( name, name, name )
            else:

                # pointer access possible for single dimensional arrays 
                # use the pointer arithmetics for accessing the pdfs
                code_str += '\tdouble %s = *( %s_flat + i_%s );\n' % ( name, name, name )

        if len( self.rv_keys ) > 0:
            code_str += self.C_code_qg
            code_str += self.rf.C_code + \
                       '// Store the values in the grid\n' + \
                       '\tmu_q +=  q * pdf;\n'
        else:
            code_str += self.rf.rf_code_eps_loop + \
                       '\tmu_q += q;\n'

        # close the random loops
        #
        for name in self.rv_keys:
            code_str += '};\n'

        if self.compiled_eps_loop:
            if self.cached_qg: # blitz matrix
                code_str += 'mu_q_arr(i_eps) = mu_q;\n'
            else:
                code_str += '*(mu_q_arr + i_eps) = mu_q;\n'
            code_str += '};\n'
        else:
            code_str += 'return_val = mu_q;'
        return code_str

    eps_arr = Property( depends_on='eps_change' )
    @cached_property
    def _get_eps_arr( self ):
        
        n_eps = self.n_eps
        min_eps = self.min_eps
        max_eps = self.max_eps
        
        # if the array of control variable points is not given
        if len( self.cv ) == 0:
            return linspace( min_eps, max_eps, n_eps )
        else:
            return self.cv 

    def _eval( self ):
        '''Evaluate the integral based on the configuration of algorithm.
        '''
        
        if self.cached_qg == False and self.compiled_qg_loop == False:
            raise NotImplementedError, \
                'Configuration for pure Python integration is too slow and is not implemented'

        self._set_compiler()
        # prepare the array of the control variable discretization
        #
        eps_arr = self.eps_arr
        mu_q_arr = zeros_like( eps_arr )

        # prepare the parameters for the compiled function in 
        # a separate dictionary
        c_params = {}

        if self.compiled_eps_loop:

            # for compiled eps_loop the whole input and output array must be passed to c
            #
            c_params['e_arr'] = eps_arr
            c_params['mu_q_arr'] = mu_q_arr
            #c_params['n_eps' ] = n_eps

        if self.compiled_qg_loop:

            # prepare the lengths of the arrays to set the iteration bounds
            #
            for rv in self.rv_list:
                c_params[ '%s_flat' % rv.name ] = rv.theta_arr

        if len( self.rv_list ) > 0:
            if self.cached_qg:
                c_params[ 'pdf_theta_grid' ] = self.pdf_theta_grid
            else:
                for rv in self.rv_list:
                    c_params['%s_pdf' % rv.name] = rv.pdf_theta_arr

        if self.cached_qg:
            conv = converters.blitz
        else:
            conv = converters.default

        t = time.time()

        if self.compiled_eps_loop:

            # C loop over eps, all inner loops must be compiled as well
            #
            from threading import Thread, _MainThread
            n_thread = 1
            C_code1 = '''for( int i_eps = 0; i_eps < 40; i_eps++){
    double eps = *( e_arr + i_eps );
double mu_q(0);
double q(0);
#line 100
double fu = 1.2e+18;
double f = 0.01;
for( int i_A = 0; i_A < 10; i_A++){
    double A = *( A_flat + i_A );
for( int i_E_mod = 0; i_E_mod < 10; i_E_mod++){
    double E_mod = *( E_mod_flat + i_E_mod );
for( int i_L = 0; i_L < 10; i_L++){
    double L = *( L_flat + i_L );
for( int i_phi = 0; i_phi < 10; i_phi++){
    double phi = *( phi_flat + i_phi );
for( int i_qf = 0; i_qf < 10; i_qf++){
    double qf = *( qf_flat + i_qf );
for( int i_z = 0; i_z < 10; i_z++){
    double z = *( z_flat + i_z );
    double pdf =  *( A_pdf + i_A)* *( E_mod_pdf + i_E_mod)* *( L_pdf + i_L)* *( phi_pdf + i_phi)* *( qf_pdf + i_qf)* *( z_pdf + i_z);

            double w = eps;
            double Le = L / 2. - z / cos( phi );
            double w_deb = exp( f * phi ) * qf * pow(Le,2.0) / E_mod / A;
            double P_deb_full = sqrt( w * E_mod * A * qf ) * exp( f * phi );
            double P_deb;
            
            // Heaviside
            if ( Le < 0 || P_deb_full > fu * A || w > w_deb ){
                P_deb = 0;
            }else{
                P_deb =P_deb_full;
            }
            
            double P_pull_x = ( Le * qf - Le * qf / ( Le - w_deb ) * ( w - w_deb ) ) * exp( f * phi );
            double P_pull;
            
            // Heaviside 
            if ( P_pull_x < 0 || w_deb > w ){
                P_pull = 0;
            }else{
                P_pull = P_pull_x;
            }
            
            // Computation of the q( ... ) function
            q = P_deb + P_pull;
        // Store the values in the grid
    mu_q +=  q * pdf;
};
};
};
};
};
};
*(mu_q_arr + i_eps) = mu_q;
};'''
            C_code2 = '''for( int i_eps = 40; i_eps < 80; i_eps++){
    double eps = *( e_arr + i_eps );
double mu_q(0);
double q(0);
#line 100
double fu = 1.2e+18;
double f = 0.01;
for( int i_A = 0; i_A < 10; i_A++){
    double A = *( A_flat + i_A );
for( int i_E_mod = 0; i_E_mod < 10; i_E_mod++){
    double E_mod = *( E_mod_flat + i_E_mod );
for( int i_L = 0; i_L < 10; i_L++){
    double L = *( L_flat + i_L );
for( int i_phi = 0; i_phi < 10; i_phi++){
    double phi = *( phi_flat + i_phi );
for( int i_qf = 0; i_qf < 10; i_qf++){
    double qf = *( qf_flat + i_qf );
for( int i_z = 0; i_z < 10; i_z++){
    double z = *( z_flat + i_z );
    double pdf =  *( A_pdf + i_A)* *( E_mod_pdf + i_E_mod)* *( L_pdf + i_L)* *( phi_pdf + i_phi)* *( qf_pdf + i_qf)* *( z_pdf + i_z);

            double w = eps;
            double Le = L / 2. - z / cos( phi );
            double w_deb = exp( f * phi ) * qf * pow(Le,2.0) / E_mod / A;
            double P_deb_full = sqrt( w * E_mod * A * qf ) * exp( f * phi );
            double P_deb;
            
            // Heaviside
            if ( Le < 0 || P_deb_full > fu * A || w > w_deb ){
                P_deb = 0;
            }else{
                P_deb =P_deb_full;
            }
            
            double P_pull_x = ( Le * qf - Le * qf / ( Le - w_deb ) * ( w - w_deb ) ) * exp( f * phi );
            double P_pull;
            
            // Heaviside 
            if ( P_pull_x < 0 || w_deb > w ){
                P_pull = 0;
            }else{
                P_pull = P_pull_x;
            }
            
            // Computation of the q( ... ) function
            q = P_deb + P_pull;
        // Store the values in the grid
    mu_q +=  q * pdf;
};
};
};
};
};
};
*(mu_q_arr + i_eps) = mu_q;
};'''

            #Thread( target=inline, args=( self.C_code, self.arg_list ), kwargs={'local_dict':c_params,
            #                              'type_converters':conv, 'verbose':0, 'compiler':'mingw32' } ).start()
#            Thread( target=inline, args=( C_code1, self.arg_list ), kwargs={'local_dict':c_params,
#                                          'type_converters':conv, 'verbose':0, 'compiler':'mingw32' } ).start()
#            print 'start1'
#            Thread( target=inline, args=( C_code2, self.arg_list ), kwargs={'local_dict':c_params,
#                                          'type_converters':conv, 'verbose':0, 'compiler':'mingw32' } ).start()
#            print 'start2'
            from multiprocessing import Process
            Process( target=inline, args=( C_code1, self.arg_list ), kwargs={'local_dict':c_params,
                                          'type_converters':conv, 'verbose':0, 'compiler':'mingw32' } ).start()
            print 'start1'
            
            Process( target=inline, args=( C_code2, self.arg_list ), kwargs={'local_dict':c_params,
                                          'type_converters':conv, 'verbose':0, 'compiler':'mingw32' } ).start()
            print 'start2'
            #inline( self.C_code, self.arg_list, local_dict=c_params, type_converters=conv, verbose=0, compiler='mingw32' ) 

        else:

            # Python loop over eps
            #
            for idx, e in enumerate( eps_arr ):
                if self.compiled_qg_loop:

                    # C loop over random dimensions
                    #
                    c_params['e'] = e # prepare the parameter
                    mu_q = inline( self.C_code, self.arg_list, local_dict=c_params,
                                   type_converters=conv, verbose=0, compiler='mingw32' )
                else:
                    # Numpy loops over random dimensions
                    #
                    # get the rf grid for all combinations of
                    # parameter values
                    #         
                    q_grid = self.rf( e, **self.param_dict )

                    # multiply the response grid with the contributions
                    # of pdf distributions (weighted by the delta of the
                    # random variable disretization)
                    #
                    q_grid *= self.pdf_theta_grid

                    # sum all the values to get the integral 
                    mu_q = sum( q_grid )

                # add the value to the return array
                mu_q_arr[idx] = mu_q

        duration = time.time() - t

        return  mu_q_arr, duration

    def _eval_mu_q( self ):
        # configure eval and call it
        pass

    def _eval_stdev_q( self ):
        # configure eval and call it
        pass

    #--------------------------------------------------------------------------------------------
    # Numpy implementation
    #--------------------------------------------------------------------------------------------
    def get_rf( self, eps ):
        '''
        Numpy based evaluation of the response function.
        '''
        return self.rf( eps, **self.param_dict )

    #---------------------------------------------------------------------------------------------
    # Output properties
    #---------------------------------------------------------------------------------------------

    # container for the data obtained in the integration
    #
    # This is not only the mean curve but also stdev and 
    # execution statistics. Such an implementation 
    # concentrates the critical part of the algorithmic 
    # evaluation and avoids duplication of code and 
    # repeated calls. The results are cached in the tuple.
    # They are accessed by the convenience properties defined
    # below. 
    #  
    results = Property( depends_on='rf_change, rand_change, conf_change, eps_change, +eps_range' )
    @cached_property
    def _get_results( self ):
        return self._eval()

    #---------------------------------------------------------------------------------------------
    # Output accessors
    #---------------------------------------------------------------------------------------------
    # the properties that access the cached results and give them a name

    mu_q_arr = Property()
    def _get_mu_q_arr( self ):
        return self.results[0]

    exec_time = Property()
    def _get_exec_time( self ):
        '''Execution time of the last evaluation.
        '''
        return self.results[1]

    mean_curve = Property()
    def _get_mean_curve( self ):
        '''Mean response curve.
        '''
        return MFnLineArray( xdata=self.eps_arr, ydata=self.mu_q_arr )

    mu_q_peak_idx = Property()
    def _get_mu_q_peak_idx( self ):
        '''Get mean peak response value'''
        return argmax( self.mu_q_arr )

    mu_q_peak = Property()
    def _get_mu_q_peak( self ):
        '''Get mean peak response value'''
        return self.mu_q_arr[ self.mu_q_peak_idx ]

    eps_at_peak = Property()
    def _get_eps_at_peak( self ):
        '''Get strain at maximum middle response mu_q
        '''
        return self.eps_arr[ self.mu_q_peak_idx ]

    stdev_mu_q_peak = Property()
    def _get_stdev_mu_q_peak( self ):
        '''
        Numpy based evaluation of the time integral.
        '''
        mu_q_peak = self.mu_q_peak
        eps_at_peak = self.eps_at_peak

        q_quad_grid = self.get_rf( eps_at_peak ) ** 2
        q_quad_grid *= self.pdf_theta_grid
        q_quad_peak = sum( q_quad_grid )
        stdev_mu_q_peak = sqrt( q_quad_peak - mu_q_peak ** 2 )

        return stdev_mu_q_peak

    #---------------------------------------------------------------------------------------------
    # Auxiliary methods
    #---------------------------------------------------------------------------------------------
    def _set_compiler( self ):
        '''Catch eventual mismatch between scipy.weave and compiler 
        '''
        try:
            uname = os.uname()[3]
        except:
            # it is not Linux - just let it go and suffer
            return

        operating_system = split( uname )[0]

        # handle the bug in the earlier versions of scipy
        # on ubuntu karmic that does not compile the 
        # converters between blitz and numpy arrays.
        #
        if operating_system == '#36-Ubuntu':
            os.environ['CC'] = 'gcc-4.1'
            os.environ['CXX'] = 'g++-4.1'


if __name__ == '__main__':

    from matplotlib import pyplot as plt

    s = SPIRRID( rf=Filament(), max_eps=0.05, n_eps=80 )

    s.add_rv( 'xi', distribution='weibull_min', scale=0.02, shape=10., n_int=30 )
    s.add_rv( 'theta', distribution='uniform', loc=0.0, scale=0.01, n_int=30 )

    s.compiled_eps_loop = False
    s.cached_qg = False

    s.mean_curve.plot( plt, color='b' , linewidth=3 )

    print s.eps_at_peak
    print s.mu_q_peak
    print s.stdev_mu_q_peak

    plt.errorbar( s.eps_at_peak, s.mu_q_peak, s.stdev_mu_q_peak, color='r', linewidth=2 )

    plt.show()
