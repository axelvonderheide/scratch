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

from enthought.traits.api import \
    HasTraits, Instance, List, Property, Array, Int, Any, cached_property, Dict, \
    Event, on_trait_change, Bool, Float

from numpy import \
    ogrid, trapz, linspace, array, sum, arange, zeros_like, \
    ones_like, max, argmax, sqrt

from i_rf import \
    IRF

from rf_filament import \
    Filament

from scipy.weave import \
    inline, converters

from stats.pdistrib.pdistrib import \
    PDistrib

from string import \
    split

import os

import time

from mathkit.mfn.mfn_line.mfn_line import \
    MFnLineArray

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
    rand_params = Dict

    # add a random variable
    # @todo change the name to add_random 
    #       and define also del_random
    def random( self, variable, distribution = 'uniform',
                loc = 0., scale = 1., shape = 1., n_int = 30 ):
        '''Declare a variable as random 
        '''
        # @todo - define the class RVariable - to avoid tuples
        #
        # @todo - check the existence of the RVariable in the params
        # 
        rand_var = PDistrib( distr_choice = distribution )
        rand_var.distr_type.set( scale = scale, shape = shape, loc = loc )
        self.rand_params[variable] = ( rand_var, n_int )

    #---------------------------------------------------------------------------------------------
    # Range of the control process variable epsilon
    #---------------------------------------------------------------------------------------------

    min_eps = Float( 0.0, eps_range = True )
    max_eps = Float( 0.05, eps_range = True )
    n_eps = Float( 80, eps_range = True )

    #--------------------------------------------------------------------
    # Define which changes in the response function and in the 
    # statistical parameters are relevant for reevaluation of the response
    #--------------------------------------------------------------------
    rf_change = Event
    @on_trait_change( 'rf.+distr' )
    def _set_rf_change( self ):
        self.rf_change = True

    rand_change = Event
    @on_trait_change( 'rand_params' )
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
    default_params = Property( Dict, depends_on = 'rf_change' )
    def _get_default_params( self ):
        '''Gather all the traits with the metadata distr specified.
        '''
        traits = self.rf.traits( distr = lambda x: x != None )
        param_dict = {}
        for tname in traits:
            param_dict[tname] = getattr( self.rf, tname )
        return param_dict

    # Discretized statistical domain
    # 
    theta_ogrid = Property( depends_on = 'rf_change, rand_change' )
    @cached_property
    def _get_theta_ogrid( self ):

        rp = self.rand_params
        ogrid_slices = [ slice( rp[pd][0].range[0], rp[pd][0].range[1], complex( 0, rp[pd][1] ) )
                     for pd in rp ]

        return ogrid[ tuple( ogrid_slices ) ]

    #---------------------------------------------------------------------------------
    # Parameters with the discretized random domain
    #---------------------------------------------------------------------------------
    params = Property( Dict, depends_on = 'rf_change, rand_change' )
    @cached_property
    def _get_params( self ):
        params = {}
        for key, val in self.default_params.items():
            params[ key ] = val

        rp = self.rand_params
        for i, name in enumerate( rp.keys() ):
            params[ name ] = self.theta_ogrid[i]
        return params

    #---------------------------------------------------------------------------------
    # PDF arrays oriented in enumerated dimensions - broadcasting possible
    #---------------------------------------------------------------------------------
    pdf_ogrid = Property( depends_on = 'rf_change, rand_change' )
    @cached_property
    def _get_pdf_ogrid( self ):
        rp = self.rand_params
        theta_ogrid = self.theta_ogrid
        return  [ rp[pd][0].distr_type.distr.pdf( theta )
                  for pd, theta in zip( rp, theta_ogrid ) ]

    #---------------------------------------------------------------------------------
    # PDF grid - mutually multiplied arrays of PDF
    #---------------------------------------------------------------------------------
    pdf_theta_grid = Property( depends_on = 'rf_change, rand_change' )
    @cached_property
    def _get_pdf_theta_grid( self ):
        return reduce( lambda x, y: x * y, self.pdf_theta_ogrid )

    #---------------------------------------------------------------------------------
    # PDF * Theta arrays oriented in enumerated dimensions - broadcasting possible
    #---------------------------------------------------------------------------------
    pdf_theta_ogrid = Property( depends_on = 'rf_change, rand_change' )
    @cached_property
    def _get_pdf_theta_ogrid( self ):
        rp = self.rand_params
        theta_ogrid = self.theta_ogrid
        pdf_theta_ogrid = [ rp[pd][0].distr_type.distr.pdf( theta ) *
                     ( theta.flatten()[1] - theta.flatten()[0] )
                     for pd, theta in zip( rp, theta_ogrid ) ]
        # cutting of the half of the first and last weight
        for arr in pdf_theta_ogrid:
            arr.ravel()[0] *= 0.5
            arr.ravel()[-1] *= 0.5

        return pdf_theta_ogrid

    #---------------------------------------------------------------------------------
    # PDF * Theta dictionary
    #---------------------------------------------------------------------------------
    # discrete values of the cdf 
    pdf_theta_dict = Property( depends_on = 'rf_change, rand_change' )
    @cached_property
    def _get_pdf_theta_dict( self ):

        theta_ogrid = self.theta_ogrid
        rp = self.rand_params

        # Memory independent = MI
        pdf_theta_dict = {}
        for pd, theta in zip( rp, theta_ogrid ):
            pdf_theta_dict[pd] = ( rp[pd][0].distr_type.distr.pdf( theta ) *
                            ( theta.flatten()[1] - theta.flatten()[0] ) ).flatten()

        # cutting of the half of the first and last weight       
        for pd in pdf_theta_dict:
            pdf_theta_dict[pd][ 0] *= .5
            pdf_theta_dict[pd][-1] *= .5

        return pdf_theta_dict

    #------------------------------------------------------------------------------------
    # Configuration of the algorithm
    #------------------------------------------------------------------------------------
    # 
    # cached_pdf_theta_grid:
    # If set to True, the cross product between the pdf values of all random variables
    # will be precalculated and stored in an n-dimensional grid
    # otherwise the product is performed for every epsilon in the inner loop anew
    # 
    cached_dG = Bool( True, alg_option = True )

    # compiled_eps_loop:
    # If set True, the loop over the control variable epsilon is compiled
    # otherwise, python loop is used.
    compiled_eps_loop = Bool( True, alg_option = True )

    # compiled_QdG_loop:
    # If set True, the integration loop over the product between the response function
    # and the pdf . theta product is performed in c
    # otherwise the numpy arrays are used.
    compiled_QdG_loop = Bool( True, alg_option = True )

    arg_list = Property( depends_on = 'rf_change, rand_change, conf_change' )
    @cached_property
    def _get_arg_list( self ):

        arg_list = []
        # create argument string for inline function
        if self.compiled_eps_loop:
            arg_list += [ 'mu_q_arr', 'e_arr' ]
        else:
            arg_list.append( 'e' )

        for name in self.params.keys():
            arg_list += ['%s_flat' % name ]

        if self.cached_dG:
            arg_list += [ 'pdf_theta_grid' ]
        else:
            arg_list += [ '%s_pdf' % name for name in self.rand_params ]

        return arg_list

    C_code_qg = Property( depends_on = 'rf_change, rand_change, conf_change' )
    @cached_property
    def _get_C_code_qg( self ):
        if self.cached_dG: # q_g - blitz matrix used to store the grid
            code_str = '\tdouble pdf = pdf_theta_grid(' + \
                       ','.join( [ 'i_%s' % name
                                  for name in self.rand_params ] ) + \
                       ');\n'
        else: # qg
            code_str = '\tdouble pdf = ' + \
                       '*'.join( [ ' *( %s_pdf + i_%s)' % ( name, name )
                                  for name in self.rand_params ] ) + \
                       ';\n'
        return code_str

    #------------------------------------------------------------------------------------
    # Configurable generation of C-code for mean curve evaluation
    #------------------------------------------------------------------------------------
    C_code = Property( depends_on = 'rf_change, rand_change, conf_change, eps_change' )
    @cached_property
    def _get_C_code( self ):

        code_str = ''
        if self.compiled_eps_loop:

            # create code string for inline function
            #
            code_str += 'for( int i_eps = 0; i_eps < %g; i_eps++){\n' % self.n_eps

            if self.cached_dG:

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

        for name, param in self.params.items():

            p = array( param ).flatten()
            n_int = len( p )
            # create the loop over the random variable
            #
            if n_int == 1:
                code_str += 'double %s = %g;\n' % ( name, param )
            else:
                code_str += 'for( int i_%s = 0; i_%s < %g; i_%s++){\n' % ( name, name, n_int, name )
                if self.cached_dG:

                    # multidimensional index needed for pdf_grid - use blitz arrays
                    #
                    code_str += '\tdouble %s = %s_flat( i_%s );\n' % ( name, name, name )
                else:

                    # pointer access possible for single dimensional arrays 
                    # use the pointer arithmetics for accessing the pdfs
                    code_str += '\tdouble %s = *( %s_flat + i_%s );\n' % ( name, name, name )

        if len( self.rand_params ) > 0:
            code_str += self.C_code_qg
            code_str += self.rf.rf_code_eps_loop + \
                       '// Store the values in the grid\n' + \
                       '\tmu_q +=  q * pdf;\n'
        else:
            code_str += self.rf.rf_code_eps_loop + \
                       '\tmu_q += q;\n'

        # close the random loops
        #
        for name, param in self.params.items():
            p = array( param ).flatten()
            n_int = len( p )
            if n_int > 1:
                code_str += '};\n'

        if self.compiled_eps_loop:
            if self.cached_dG: # blitz matrix
                code_str += 'mu_q_arr(i_eps) = mu_q;\n'
            else:
                code_str += '*(mu_q_arr + i_eps) = mu_q;\n'
            code_str += '};\n'
        else:
            code_str += 'return_val = mu_q;'
        return code_str

    def _eval( self ):
        '''Evaluate the integral based on the configuration of algorithm.
        '''
        print 'REEVALUATION'

        if self.cached_dG == False and self.compiled_QdG_loop == False:
            raise NotImplementedError, \
                'Configuration for pure Python integration is too slow and is not implemented'

        n_eps = self.n_eps
        min_eps = self.min_eps
        max_eps = self.max_eps

        self._set_compiler()
        # prepare the array of the control variable discretization
        #
        eps_arr = linspace( min_eps, max_eps, n_eps )
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

        if self.compiled_QdG_loop:

            # prepare the lengths of the arrays to set the iteration bounds
            #
            for name, param in self.params.items():
                p = array( param ).flatten()
                #c_params[ 'n_%s' % name ] = len( p )
                c_params[ '%s_flat' % name ] = p

        if len( self.rand_params ) > 0:
            if self.cached_dG:
                c_params[ 'pdf_theta_grid' ] = self.pdf_theta_grid
            else:
                for name in self.rand_params:
                    c_params['%s_pdf' % name] = self.pdf_theta_dict[name]

        if self.cached_dG:
            conv = converters.blitz
        else:
            conv = converters.default

        t = time.clock()

        if self.compiled_eps_loop:

            # C loop over eps, all inner loops must be compiled as well
            #
            inline( self.C_code, self.arg_list, local_dict = c_params,
                    type_converters = conv, verbose = 0 )

        else:

            # Python loop over eps
            #
            for idx, e in enumerate( eps_arr ):

                if self.compiled_QdG_loop:

                    # C loop over random dimensions
                    #
                    c_params['e'] = e # prepare the parameter
                    mu_q = inline( self.C_code, self.arg_list, local_dict = c_params,
                                   type_converters = conv, verbose = 0 )
                else:

                    # Numpy loops over random dimensions
                    #
                    # get the rf grid for all combinations of
                    # parameter values
                    #                    
                    q_grid = self.rf( e, **self.params )

                    # multiply the response grid with the contributions
                    # of pdf distributions (weighted by the delta of the
                    # random variable disretization)
                    #
                    q_grid *= self.pdf_theta_grid

                    # sum all the values to get the integral 
                    mu_q = sum( q_grid )

                # add the value to the return array
                mu_q_arr[idx] = mu_q

        duration = time.clock() - t

        return  eps_arr, mu_q_arr, duration

    #--------------------------------------------------------------------------------------------
    # Numpy implementation
    #--------------------------------------------------------------------------------------------
    def get_pq( self, eps ):
        '''
        Numpy based evaluation of the response function.
        '''
        print 'calling pq with', self.params
        return self.rf( eps, **self.params )

    def eval_peqg( self, min_eps, max_eps, n_eps ):
        '''
        Numpy based evaluation of the time integral.
        '''
        eps_arr = linspace( min_eps, max_eps, n_eps )

        mu_q = zeros_like( eps_arr )
        t = time.clock()
        for id, eps in enumerate( eps_arr ):
            q_grid = self.get_pq( eps )
            # why cdf?
            for cdf in self.pdf_theta_ogrid:
                q_grid *= cdf
            mu_q[id] = sum( q_grid )

        duration = time.clock() - t

        # standard deviation at peak force

        mu_max = max( mu_q )
        eps_max = eps_arr[argmax( mu_q )]
        print 'eps_max', eps_max
        q_quad_grid = self.get_pq( eps_max ) ** 2
        print 'q_quad_max', sum( q_quad_grid )
        for weighted_pdf in self.pdf_theta_ogrid:
            q_quad_grid *= weighted_pdf
        q_quad_max = sum( q_quad_grid )
        print 'q_quad_max', q_quad_max
        print 'mu_max', mu_max ** 2
        stdev = sqrt( q_quad_max - mu_max ** 2 )

        return eps_arr, mu_q, duration, eps_max, stdev

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
    results = Property( depends_on = 'rf_change, rand_change, conf_change, eps_change' )
    @cached_property
    def _get_results( self ):
        return self._eval()

    #---------------------------------------------------------------------------------------------
    # Output accessors
    #---------------------------------------------------------------------------------------------
    # the properties that access the cached results and give them a name
    mean_curve = Property()
    @cached_property
    def _get_mean_curve( self ):
        '''Mean response curve.
        '''
        return MFnLineArray( xdata = self.results[0], ydata = self.results[1] )

    exec_time = Property()
    def _get_exec_time( self ):
        '''Execution time of the last evaluation.
        '''
        return self.results[2]

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
    s = SPIRRID( rf = Filament(), max_eps = 0.05, n_eps = 80 )


#    s.random( 'xi', distribution = 'weibull_min', scale = 0.02, shape = 10. )
#    print 'params', s.params
#    print s.C_code
#    eps, mu_q, d = s.eval( 0.0, 0.05, 20 )
#    print 'mu_q', mu_q
#    print 'duration', d

#    code_str, arg_list = s.get_ceqg_pointer()
#    print code_str    
    #X, Y = s.eval()
    #plt.plot( X, Y, linewidth=2, color='red' )
    #plt.show()

    from matplotlib import pyplot as plt

    s.random( 'theta', distribution = 'uniform', loc = 0.0, scale = 0.01, n_int = 30 )

    print s.C_code
    print s.pdf_theta_grid

#    s.mean_curve.plot( plt, color = 'r' , linewidth = 3 )

    X, Y, dur, eps_max, stdev = s.eval_peqg( min_eps = 0.0, max_eps = .05, n_eps = 80 )
    plt.plot( X, Y, linewidth = 2, color = 'b' )
    plt.errorbar( eps_max, max( Y ), stdev, color = 'r', linewidth = 2 )
    plt.show()
