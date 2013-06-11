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
# Created on Apr 9, 2010 by: rch

from enthought.traits.api import \
    HasTraits, Array, Property, cached_property

import pylab as p

from stats.pdistrib.pdistrib import \
    IPDistrib, PDistrib

from numpy import \
    ogrid, frompyfunc, array, \
    linspace, hstack, vstack, vectorize, \
    apply_over_axes, sum, where, ones_like, \
    zeros_like, arange, zeros, sign

from numpy.random import \
    random_integers

from scipy.integrate import \
    trapz, simps, romb

from scipy import \
    weave

from math import pi

import os

import time

# Quantities for the response function
# and randomization
# 
E_mod = 70 * 1e+9 # Pa
sig_u = 1.25 * 1e+9 # Pa
D = 26 * 1.0e-6 # m
A = ( D / 2.0 ) ** 2 * pi
xi_u = sig_u / E_mod

os.environ['CC'] = 'gcc-4.1 -O2'
os.environ['CXX'] = 'g++-4.1 -O2'

#---------------------------------------------------------------------------
# RF - VERSION A - SCALAR RESPONSE FUNCTION WITH SUBSEQUENT VECTORIZATION
#---------------------------------------------------------------------------
def qA_fn( eps, xi, theta, lambda_, A ):
    '''
    evaluate the response function for all combinations of [eps,xi,theta] 
    (see Part 1, Stochastic Modeling of Multifilament Yarns, 06)
    '''
    # eps is global and eps_ is local strain
    eps_ = ( ( eps - theta * ( 1 + lambda_ ) ) /
            ( ( 1 + theta ) * ( 1 + lambda_ ) ) )
    if eps_ < 0 or eps_ > xi :
        return 0.0
    else:
        EA = E_mod * A
        return EA * eps_

# vectorized version
qA_vct_fn = vectorize( qA_fn, [float] )

def Heaviside( x ):
    return ( sign( x ) + 1.0 ) / 2.0

    #return sign( sign( x ) + 1 )


#--------------------------------------------------------------------------
# RF - VERSION B: RESPONSE FUNCTION USING ARRAY EXPRESSIONS
#--------------------------------------------------------------------------
def qB_vct_fn( eps, xi, theta, lambda_, A ):
    '''
    Implements the response function with arrays as variables.
    first extract the variable discretizations from the orthogonal grid.
    '''

    # NOTE: as each variable is an array oriented in different direction
    # the algebraic expressions (-+*/) perform broadcasting,. i.e. performing
    # the operation for all combinations of values. Thus, the resulgin eps
    # is contains the value of local strain for any combination of 
    # global strain, xi, theta and lambda 
    #
    eps_ = ( eps - theta * ( 1 + lambda_ ) ) / ( ( 1 + theta ) * ( 1 + lambda_ ) )

    # cut off all the negative strains due to delayed activation
    # 
    eps_ *= Heaviside( eps_ )

    # broadcast eps also in the xi - dimension 
    # (by multiplying with array containing ones with the same shape as xi )
    #
    eps_grid = eps_ * Heaviside( xi - eps_ )

    # cut off all the realizations with strain greater than the critical one.
    # 
    # eps_grid[ eps_grid >= xi ] = 0

    # transform it to the force
    # 
    q_grid = E_mod * A * eps_grid

    return q_grid

#---------------------------------------------------------------------------
# RF - VERSION C - SCALAR RESPONSE FUNCTION WITH SUBSEQUENT VECTORIZATION
# using WEAVE compilation into C++ code.
#---------------------------------------------------------------------------
# vectorize the response function to work on arrays
# it gets three arrays as input (strain, xi, theta, lambda, A)
#
# REMARK: the slowest version - it does not bring anything
# due to the overhead of calling the compiled function is
# too large.
#
code = '''
    double eps_;
    double q;
    double xlambda( lambd );
    double xtheta( theta );
    double xeps( eps );
    double xA( A );
    double xxi( xi );
    eps_ = ( ( xeps - xtheta * ( 1.0 + xlambda ) ) / 
            ( ( 1.0 + xtheta ) * ( 1.0 + xlambda ) ) );
    if ( eps_ < 0 || eps_ > xxi ){
        q = 0.0;
    }else{
        q = %g * xA * eps_;
    }
    return_val = q;
    ''' % E_mod

def qC_weave_fn( eps, xi, theta, lambd, A ):
    return weave.inline( code, ['eps', 'xi', 'theta', 'lambd', 'A'], verbose = 0 )

qC_vct_fn = vectorize( qC_weave_fn, [float] )

#--------------------------------------------------------------------------
# RF - VERSION D: RESPONSE FUNCTION USING SLICED EXPRESSIONS
#--------------------------------------------------------------------------

class qD_vct_fn( HasTraits ):
    '''Numpy version avoiding copying of interim results
    '''
    pdf_grid = Array

    def __init__( self, pdf_grid, **kw ):
        super( qD_vct_fn, self ).__init__( **kw )
        self.pdf_grid = pdf_grid

    eps_grid = Property( Array, depends_on = 'pdf_grid' )
    @cached_property
    def _get_eps_grid( self ):
        return zeros_like( self.pdf_grid )

    q_grid = Property( Array, depends_on = 'pdf_grid' )
    @cached_property
    def _get_q_grid( self ):
        return zeros_like( self.pdf_grid )

    def __call__( self, eps, xi, theta, lambda_, A ):
        '''
        @todo: implement the sliced array expressions. This avoids
        the allocation of memory for interim results and can be 
        "weaved and blitzed" in a reasonable way.
        '''

        # NOTE: as each variable is an array oriented in different direction
        # the algebraic expressions (-+*/) perform broadcasting,. i.e. performing
        # the operation for all combinations of values. Thus, the resulgin eps
        # is contains the value of local strain for any combination of 
        # global strain, xi, theta and lambda 
        #
        self.eps_grid[...] = ( ( eps - theta *
                              ( 1 + lambda_ ) ) /
                              ( ( 1 + theta ) * ( 1 + lambda_ ) ) )

        self.eps_grid[ self.eps_grid < 0 ] = 0
        self.eps_grid[ self.eps_grid >= xi ] = 0

        self.q_grid[...] = E_mod * A * self.eps_grid

        return self.q_grid

#---------------------------------------------------------------------------
# RF - VERSION E - VECTOR RESPONSE FUNCTION WITH LOOPS COMPILED
# using WEAVE compilation into C++ code.
#---------------------------------------------------------------------------
#
class qE_vct_fn( HasTraits ):
    '''Numpy version avoiding copying of interim results
    '''
    pdf_grid = Array

    def __init__( self, pdf_grid, **kw ):
        super( qE_vct_fn, self ).__init__( **kw )
        self.pdf_grid = pdf_grid

    def integ( self, eps, Xi, Theta, Lambda, Area ):

        Xi_flat = Xi.flatten()
        Theta_flat = Theta.flatten()
        Lambda_flat = Lambda.flatten()
        Area_flat = Area.flatten()

        n_xi = len( Xi_flat )
        n_theta = len( Theta_flat )
        n_lambda = len( Lambda_flat )
        n_A = len( Area_flat )

        pdf_grid = self.pdf_grid

        code = '''
        double xeps = eps;
        double eps_;
        double q;

        double d_xi = Xi_flat( 1 ) - Xi_flat( 0 ); 
        double d_theta = Theta_flat( 1 ) - Theta_flat( 0 );  
        double d_lambda = Lambda_flat( 1 ) - Lambda_flat( 0 );  
        double d_A = Area_flat( 1 ) - Area_flat( 0 );
        
        double d_Theta = d_xi * d_theta * d_lambda * d_A;

        double mu_q = 0;
        
        for( int i_xi = 0; i_xi < n_xi; i_xi++ ){
        
            double xi = Xi_flat( i_xi );
            
            for( int i_theta = 0; i_theta < n_theta; i_theta++ ){
            
                double theta = Theta_flat( i_theta );
                
                for( int i_lambda = 0; i_lambda < n_lambda; i_lambda++ ){
                
                    double lambda = Lambda_flat( i_lambda );
                    
                    for( int i_A = 0; i_A < n_A; i_A++ ){
                    
                        double A = Area_flat( i_A );
                    
                        // Computation of the q( ... ) function
                        //
                        eps_ = ( ( xeps - theta * ( 1.0 + lambda ) ) / 
                                ( ( 1.0 + theta ) * ( 1.0 + lambda ) ) );
                        if ( eps_ < 0 || eps_ > xi ){
                            q = 0.0;
                        }else{
                            q = %g * A * eps_;
                        }
                        
                        // Get the current pdf
                        //
                        double pdf = pdf_grid( i_xi, i_theta, i_lambda, i_A );
                        
                        // Store the values in the grid
                        //
                        mu_q += q * pdf * d_Theta;
                    };
                };
            };
        };
        return_val = mu_q;
        ''' % E_mod

        mu_q = weave.inline( code, ['eps',
                                   'Xi_flat', 'n_xi',
                                   'Theta_flat', 'n_theta',
                                   'Lambda_flat', 'n_lambda',
                                   'Area_flat', 'n_A',
                                   'pdf_grid'
                                   ],
                            type_converters = weave.converters.blitz,
                            verbose = 0 )
        return mu_q

#--------------------------------------------------------------------------
# INTEGRATION
#--------------------------------------------------------------------------
def integ_qtheta( q_grid, pdf_grid, theta_ogrid ):
    '''Integrate the multidimensional integral 
    int( q(eps,theta) * pdf(theta), d_theta )
    '''
    # multiply the response function with the pdf grid 
    # q[eps,xi,theta,lambda,A] * pdf[xi,theta,lambda,A]
    # NOTE: pdf has no eps dimension - thus, it will be broadcasted
    # to all values of eps, i.e. the same pdf_grid will be used for
    # all strains
    #
    qpdf_grid = q_grid * pdf_grid

    # perform the integration over theta and xi and lambda_
    # NOTE: after the inner call to trapz, the resulting gets
    # reduced by one dimension, therefore, integration starts from
    # the index 3 and goes to 1
    #

    def integ_qtheta_recur( q_x_grid, theta_ogrid ):
        '''Recursive integration.
        all dimensions of theta_ogrid are integrated.
        the q_x_grid is assumed to have AT LEAST the same n_dims as theta_ogrid.
        Integration starts with the last axis. The result has the dimension
        q_x_grid.ndims - theta_ogrid.ndims
        '''
        # take the last dimension of the array 
        q_axis = q_x_grid.ndim - 1
        theta_axis = len( theta_ogrid ) - 1
        if theta_axis >= 0:
            # reduces the dimension of the array expression by one.
            q_x_grid = trapz( 
                          q_x_grid,
                          theta_ogrid[ theta_axis ].flatten(), axis = q_axis
                          )
            return integ_qtheta_recur( q_x_grid, theta_ogrid[:-1] )
        else:
            return q_x_grid

    mu_q_arr = integ_qtheta_recur( qpdf_grid, theta_ogrid )

    return mu_q_arr

#--------------------------------------------------------------------------
# PLOTTING
#--------------------------------------------------------------------------
def plot_scattered_realizations( q_vct_fn, eps_arr, theta_ogrid ):
    '''Plot only a fraction of realizations.
    '''
    # go through the arrays and extract n_xp in each theta direction
    # 
    ix_list = [ arange( 0, ix.flatten().size, ix.flatten().size / 10 ) for ix in theta_ogrid ]
    # construct the reduced space of the random variables
    #
    theta_minmax = [ th.flatten()[ix] for ix, th in zip( ix_list, theta_ogrid ) ]
    theta_minmax_ogrid = ( theta_minmax[0][:, None, None, None],
                          theta_minmax[1][None, :, None, None],
                          theta_minmax[2][None, None, :, None],
                          theta_minmax[3][None, None, None, :] )

    # and let the standard plot_realizations do the job
    #
    plot_realizations( q_vct_fn, eps_arr, theta_minmax_ogrid )

def plot_realizations( q_vct_fn, eps_arr, theta_ogrid ):
    '''Plot the realizations for all combinations of parameters
    stored in theta_ogrid.
    
    The full expansion is used here by adding the eps_arr to the
    theta_ogrid and running the q_vct_fn over the expanded grid.
    The array is then reshaped such that the q_curves can be 
    accessed by slicing the q_grid_ over the first dimension for 
    a given second dimension that corresponds to a combination
    of theta parameters.
    '''
    eps = eps_arr[:, None, None, None, None]
    eps_theta_ogrid = [ eps ] + [ theta[None, ...] for theta in theta_ogrid ]
    grid_size = reduce( lambda x, y: x * y, [ s.size for s in eps_theta_ogrid ] )
    print 'GRID SIZE: %d:' % grid_size
    if grid_size > 20000000:
        print 'That\'s too much! You have to implement it better'
        print 'or just don\'t take that many points ;-)'
        print 'you might double this when running only one version'
        return

    q_grid = q_vct_fn( *eps_theta_ogrid ) # calculate (* make tuple/list to args)

    # plot all the realizations stored in q_grid
    # first reshape the q_grid such that it has
    # epsilons in the first index and the index s1 of the 
    # random parameter combination in the second index s2.
    # (s1 and s2 was defined above)
    #
    s1 = q_grid.shape[0]
    s2 = q_grid.shape[1] * q_grid.shape[2] * q_grid.shape[3] * q_grid.shape[4]
    q_grid_ = q_grid.reshape( ( s1, s2 ) )

    if s2 > 4000:
        # this would be too much for plotting
        # take randomly some 4000 realizations
        plot_idx_list = random_integers( 0, high = s2 - 1, size = 4000 )
    else:
        # plot them all
        plot_idx_list = range( s2 )
    for s in plot_idx_list:
        q_arr = q_grid_[:, s]
        p.plot( eps_arr, q_arr, color = 'grey' )

#--------------------------------------------------------------------------
# INTEGRATION OF THE PRODUCT Q * PDF - LOOPLESS
#--------------------------------------------------------------------------
def integ_qtheta_eps_loopless( q_vct_fn, eps_arr, pdf_grid, theta_ogrid ):
    '''
    Fully expanded orthogonal grid of variables (including epsilon).

    Construct an orthogonal grid of arrays for 
    epsilon: in the first dimension (note that empty axes 1- and 2- were
        using the construct eps[:,None,None] to prepare it for
        broadcasting (implicit copy of values for any combination of
               [xi,theta]
      xi     : in the second dimension (add the first empty dimension)
               for broadcasting (the same xi is taken for any value of eps
      theta  : the same as xi - therefore the list loop can be used.
      lambda : the same as xi - therefore the list loop can be used.
    '''
    eps = eps_arr[:, None, None, None, None]
    eps_theta_ogrid = [ eps ] + [ theta[None, ...] for theta in theta_ogrid ]
    grid_size = reduce( lambda x, y: x * y, [ s.size for s in eps_theta_ogrid ] )
    print 'GRID SIZE: %d:' % grid_size
    if grid_size > 20000000:
        print 'That\'s too much! You have to implement it better'
        print 'or just don\'t take that many points ;-)'
        print 'you might double this when running only one version'
        return

    q_grid = q_vct_fn( *eps_theta_ogrid ) # calculate (* make tuple/list to args)
    mu_q_arr = integ_qtheta( q_grid, pdf_grid, theta_ogrid )

    return mu_q_arr

#--------------------------------------------------------------------------
# INTEGRATION OF THE PRODUCT Q * PDF - LOOPED OVER eps
#--------------------------------------------------------------------------
def integ_qtheta_eps_looped( q_vct_fn, eps_arr, pdf_grid, theta_ogrid ):
    '''
    Integration with python loop over the eps dimension.
    (the fastest implementation)
    '''
    mu_q_arr = zeros_like( eps_arr )
    for idx, eps in enumerate( eps_arr ):

        q_grid = q_vct_fn( eps, *theta_ogrid )
        mu_q = integ_qtheta( q_grid, pdf_grid, theta_ogrid )
        mu_q_arr[idx] = mu_q

    return mu_q_arr

#--------------------------------------------------------------------------
# INTEGRATION OF THE PRODUCT Q * PDF - LOOPED OVER eps
#--------------------------------------------------------------------------
def integ_qtheta_eps_cum( q_vct_fn, eps_arr, pdf_grid, theta_ogrid ):
    '''
    Integration with python loop over the eps dimension.
    (the fastest implementation)
    '''
    mu_q_arr = zeros_like( eps_arr )
    for idx, eps in enumerate( eps_arr ):

        mu_q = q_vct_fn.integ( eps, *theta_ogrid )
        mu_q_arr[idx] = mu_q

    return mu_q_arr


def run():

    # discretization parameters
    #
    np = 30
    n_eps = 80

    # implementation version
    plot_realizations = False
    version1A = False  # loopless with scalar python response function 
    version1B = False  # loopless with vector numpy response function
    version1C = False  # loopless with scalar weave-c response function
    version2A = False  # looped with scalar python response function 
    version2B = False  # looped with vector numpy response function
    version2C = False  # looped with scalar weave-c response function
    version2D = False  # looped with sliced numpy response function
    version2E = False   # looped with looped weave-c response function
    version2F = True   # looped with looped weave-c response function

    #-------------------------------------------------------------------
    # DEFINE THE RANDOM VARIABLES
    #-------------------------------------------------------------------
    # Corresponding to the response function define a list
    # PDistrib instances specifying the distribution and
    # its parameters.
    #
    print 'DISTR XI'
    pxi = PDistrib( distr_choice = 'weibull_min' )
    pxi.distr_type.set( scale = 0.02, shape = 10 )

    print 'DISTR THETA'
    ptheta = PDistrib( distr_choice = 'uniform' )
    ptheta.distr_type.set( loc = 0.0, scale = 0.01 )

    print 'DISTR LAMBDA'
    plambda = PDistrib( distr_choice = 'uniform' )
    plambda.distr_type.set( loc = 0.0, scale = 0.2 )

    print 'DISTR A'
    pA = PDistrib( distr_choice = 'uniform' )
    pA.distr_type.set( loc = A, scale = A * 0.3 )

    pd_list = [ pxi, ptheta, plambda, pA ]

    start = time.clock() # take the start time

    print 'PREPARING THE STATISTICAL DOMAIN'
    #-------------------------------------------------------------------
    # DISCRETIZATION OF THE RANDOM VARIABLES
    #-------------------------------------------------------------------
    # use the pdf_list to generate discretizations for each variable.
    # The discretized arrays are put into a tuple. Each array has
    # the dimension of the whole statistical domain, i.e. 4 in the
    # running example, but it has the length 1 in all dimensions 
    # different to its own dimension (or index). This job is done
    # by the numpy.ogrid generator. 
    #
    # The following code is a generic shortcut for :
    #
    #    theta_ogrid = ogrid[ pxi.range[0]:pxi.range[1]:complex( 0, np ), 
    #                         ptheta.range[0]:ptheta.range[1]:complex( 0, np ),
    #                         plambda.range[0]:plambda.range[1]:complex( 0, np ),
    #                         pA.range[0]:pA.range[1]:complex( 0, np ) ]
    #
    ogrid_slices = [ slice( pd.range[0], pd.range[1], complex( 0, np ) )
                     for pd in pd_list ]

    # by converting the list to a tuple, ogrid handles it as 
    # a normal set of comma-separated parameters as in the example above
    # 
    theta_ogrid = ogrid[ tuple( ogrid_slices ) ]

    #-------------------------------------------------------------------
    # PDFS FOR INDIVIDUAL RANDOM VARIABLE DISCRETIZATIONS
    #-------------------------------------------------------------------
    # get the pdf values for xi, theta, ... values
    # the generic code for arbitrary list of distributions is below
    # annd corresponds to the following code for the 
    # considered 4-random variables a filament response.
    # 
    #    pdf_xi     = pxi.distr_type.distr.pdf( theta_ogrid[0] )
    #    pdf_theta  = ptheta.distr_type.distr.pdf( theta_ogrid[1] )
    #    pdf_lambda = plambda.distr_type.distr.pdf( theta_ogrid[2] )
    #    pdf_A      = pA.distr_type.distr.pdf( theta_ogrid[3] )
    #
    # NOTE: the xi and theta have different orientations
    #       while pdf_xi goes along the 0-axis i.e. [:,0]
    #       the pdf_theta goes along the 1-axis i.e. [0,1]
    #       etc...

    pdf_list = [ pd.distr_type.distr.pdf( theta )
                 for pd, theta in zip( pd_list, theta_ogrid ) ]

    #-------------------------------------------------------------------
    # PDF VALUES OVER THE DISCRETIZED RANDOM DOMAIN
    #-------------------------------------------------------------------
    # multiply the discretized pdfs to get the values of pdf for
    # all combinations of xi and theta and lambda
    # NOTE: the result of multiplying [10,1] * [1,10] array is
    #       [10,10] array
    #
    # NOTE: The below code applies recursive product through
    #       the list and corresponds to the explicit version 
    #  
    # pdf_grid = pdf_xi * pdf_theta * pdf_lambda * pdf_A
    #
    pdf_grid = reduce( lambda x, y: x * y, pdf_list )

    #---------------------------------------------------------------------
    # STRAIN DISCRETIZATION
    #---------------------------------------------------------------------
    eps_arr = linspace( 0., 0.05, n_eps )

    # plot selected realizations from the random domain
    # (use qB_vct_fn, as it is faster than qA_vct_fn)
    #
    if plot_realizations:
        plot_scattered_realizations( qB_vct_fn, eps_arr, theta_ogrid )

    #----------------------------------------------------------------------------
    # TEST CALCULATION USING THE COMBINED IMPLEMENTATION OF RF AND INTEGRATION
    #----------------------------------------------------------------------------

    # test the provided immplementations
    #
    if version1A:
        print '----------------------------------------------------------------'
        print 'VERSION 1A: VECTORIZED SCALAR RESPONSE FUNCTION x LOOPLESS INTEG'
        start = time.clock() # take the start time
        mu_q_arr = integ_qtheta_eps_loopless( qA_vct_fn, eps_arr, pdf_grid, theta_ogrid )
        end = time.clock() # take the end time
        print 'Duration of response function evaluation:', end - start, 's'
        p.plot( eps_arr, mu_q_arr, color = 'yellow', linewidth = 7 )

    if version1B:
        print '---------------------------------------------------------------'
        print 'VERSION 1B: ARRAY RESPONSE FUNCTION x LOOPLESS INTEG           '
        start = time.clock() # take the start time
        mu_q_arr = integ_qtheta_eps_loopless( qB_vct_fn, eps_arr, pdf_grid, theta_ogrid )
        end = time.clock()
        print 'Duration of response function evaluation:', end - start, 's'
        p.plot( eps_arr, mu_q_arr, color = 'red', linewidth = 5 )

    if version1C:
        print '---------------------------------------------------------------'
        print 'VERSION 1C: VECTORIZED WEAVE RESPONSE FUNCTION x LOOPLESS INTEG'
        start = time.clock() # take the start time
        mu_q_arr = integ_qtheta_eps_loopless( qC_vct_fn, eps_arr, pdf_grid, theta_ogrid )
        end = time.clock()
        print 'Duration of response function evaluation:', end - start, 's'
        p.plot( eps_arr, mu_q_arr, color = 'yellow', linewidth = 5 )

    if version2A:
        print '---------------------------------------------------------------'
        print 'VERSION 2A: VECTORIZED SCALAR RESPONSE FUNCTION x LOOPED INTEG'
        start = time.clock() # take the start time
        mu_q_arr = integ_qtheta_eps_looped( qA_vct_fn, eps_arr, pdf_grid, theta_ogrid )
        end = time.clock()
        print 'Duration of response function evaluation:', end - start, 's'
        print '--------------------------------------------------------------'
        p.plot( eps_arr, mu_q_arr, color = 'blue', linewidth = 3 )

    if version2B:
        print '---------------------------------------------------------------'
        print 'VERSION 2B: ARRAY RESPONSE FUNCTION x LOOPED INTEG           '
        start = time.clock() # take the start time
        mu_q_arr = integ_qtheta_eps_looped( qB_vct_fn, eps_arr, pdf_grid, theta_ogrid )
        end = time.clock()
        print 'Duration of response function evaluation:', end - start, 's'
        print '--------------------------------------------------------------'
        p.plot( eps_arr, mu_q_arr, color = 'blue', linewidth = 3 )

    if version2C:
        print '---------------------------------------------------------------'
        print 'VERSION 2C: VECTORIZED WEAVE RESPONSE FUNCTION x LOOPED INTEG'
        start = time.clock() # take the start time
        mu_q_arr = integ_qtheta_eps_looped( qC_vct_fn, eps_arr, pdf_grid, theta_ogrid )
        end = time.clock()
        print 'Duration of response function evaluation:', end - start, 's'
        p.plot( eps_arr, mu_q_arr, color = 'yellow', linewidth = 5 )

    if version2D:
        print '---------------------------------------------------------------'
        print 'VERSION 2D: SLICED RESPONSE FUNCTION x LOOPED INTEG           '
        start = time.clock() # take the start time
        mu_q_arr = integ_qtheta_eps_looped( qD_vct_fn( pdf_grid ), eps_arr, pdf_grid, theta_ogrid )
        end = time.clock()
        print 'Duration of response function evaluation:', end - start, 's'
        print '--------------------------------------------------------------'
        p.plot( eps_arr, mu_q_arr, color = 'blue', linewidth = 3 )

    if version2E:
        print '---------------------------------------------------------------'
        print 'VERSION 2E: SLICED RESPONSE FUNCTION x LOOPED INTEG           '
        start = time.clock() # take the start time
        mu_q_arr = integ_qtheta_eps_looped( qE_vct_fn( pdf_grid ), eps_arr, pdf_grid, theta_ogrid )
        end = time.clock()
        print 'Duration of response function evaluation:', end - start, 's'
        print '--------------------------------------------------------------'
        p.plot( eps_arr, mu_q_arr, color = 'black', linewidth = 1 )

    if version2F:
        print '---------------------------------------------------------------'
        print 'VERSION 2E: SLICED RESPONSE FUNCTION x LOOPED INTEG           '
        start = time.clock() # take the start time
        print theta_ogrid[2].flatten()
        mu_q_arr = integ_qtheta_eps_cum( qE_vct_fn( pdf_grid ), eps_arr, pdf_grid, theta_ogrid )
        end = time.clock()
        print 'Duration of response function evaluation:', end - start, 's'
        print '--------------------------------------------------------------'
        p.plot( eps_arr, mu_q_arr, color = 'black', linewidth = 1 )


    # FURTHER VERSIONS
    #
    # 1)
    # Decomposing the integration domain further with looped integration of
    # over the individual random variables?
    #
    # This might be necessary for really large problems with more than 5
    # random variables and need for more than 20 points in each 
    # theta direction.
    #
    # 2)
    # scipy.weave.blitz provides a machinery how to implement loops 
    # over numpy arrays directly in c++. It might provide even more efficient
    # evaluation of RF on a compact array. However, I do not expect too 
    # much improvement compared to the array-based implementation.
    #

    # show all plots performed so far in a matplotlib window
    #
    p.show()

if __name__ == '__main__':

    #  
    use_profiling = False

    if use_profiling:
        import cProfile
        cProfile.run( 'run()', 'spirrid_integ_tprof' )

        import pstats
        p = pstats.Stats( 'spirrid_integ_tprof' )
        p.strip_dirs()
        print 'cumulative'
        p.sort_stats( 'cumulative' ).print_stats( 50 )

    else:
        run()
