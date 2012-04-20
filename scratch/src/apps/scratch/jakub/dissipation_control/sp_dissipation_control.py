'''
Created on Jun 30, 2010

@author: jakub
'''
#TODO: start load step from origin - why the dissipation control with iter 
#flag does not work?
from sympy import symbols, solve, integrate
from numpy import linspace, array, linalg
from matplotlib import pylab
from math import fabs

x, y = symbols( 'xy' )

fn = -( x - 1 ) ** 2 + 1
dfn = -2 * ( x - 1 )

tau = 5.0e-3
tol = 1.0e-6

#x1 = 0.1
#y1 = f.subs( x, x1 )

def get_val_ana( lambda_n ):
    '''
    get the x value
    @param lambda_n:
    '''
    x1 = solve( fn - lambda_n, x )[1]
    return x1

def compute_dissipation( lambda_n, lambda_k, U_n, U_k ):
    '''
    compute analytically the dissipation in arbitrary step
    lambda_n -> lambda_k
    U_n -> U_k
    @param lambda_n:
    @param lambda_k:
    @param U_n:
    @param U_k:
    '''
    U_tot = integrate( fn, ( x, 0., U_k ) )
    #print 'U_tot ', U_tot
    U_el = lambda_k * U_k / 2.
    #print 'U_el ', U_el
    U_dis = integrate( fn, ( x, 0., U_n ) ) - lambda_n * U_n / 2.
    U_inel = U_tot - U_el - U_dis
    return U_inel


def plot_fn():
    '''
    plot the fn in matplotlib
    '''
    sam = linspace( 0., 1. )#default 50
    val = []
    for i in sam:
        val.append( fn.subs( x, i ) )

    pylab.plot( sam, val )

def get_dif_d( lambda_i, U_i ):
    '''
    get derivation of dissipative fn for given lamda and displacement
    @param lambda_i:
    @param U_i:
    '''
    f_bar = 1.
    h = lambda_i * f_bar / 2.
    w = -f_bar * U_i / 2.
    return h, w

def get_dif_K( lambda_i, U_i ):
    '''
    get stiffness matrix
    @param lambda_i:
    @param U_i:
    '''
    K = dfn.subs( x, U_i )
    return K

def get_dissipation( lambda_n, lambda_k, U_n, U_k, flag ):
    '''
    
    @param lambda_n:
    @param lambda_k:
    @param U_n:
    @param U_k:
    @param flag:
    '''
    f_bar = 1.
    g = f_bar * ( lambda_n * ( U_k - U_n ) - \
              ( lambda_k - lambda_n ) * U_n )\
            / 2.
    return g


def get_val_numNR( lambda_n, lambda_k, U_n, flag = 'step' ):
    '''
    Newton-Raphson procedure 
    standart with iter flag, modified with step flag
    @param lambda_n:
    @param lambda_k:
    @param U_n:
    @param flag:
    '''
    U_k = U_n
    f_bar = 1.

    p_x = [U_n]
    p_y = [lambda_n]
    K = get_dif_K( lambda_n, U_n )

    d_U = ( lambda_k - lambda_n ) * f_bar / K
    U_k += d_U
    p_x.append( U_k )
    p_y.append( lambda_k )
    res = lambda_k * f_bar - fn.subs( x, U_k )
    p_x.append( U_k )
    p_y.append( lambda_k - res )

    while res > tol:
        if flag == 'iter':
            K = get_dif_K( lambda_k, U_k )
        d_U = res / K
        U_k += d_U
        p_x.append( U_k )
        p_y.append( lambda_k )
        res = lambda_k * f_bar - fn.subs( x, U_k )
        p_x.append( U_k )
        p_y.append( lambda_k - res )
    pylab.plot( p_x, p_y )

def get_val_numD( lambda_n, lambda_k, U_n, flag = 'step' ):
    '''
    dissipation control algorithm [Verhoosel, Eckhard]
    updating in every iteration: 
    A mtx - flag iter
    K mtx - flag K-iter
    nothing - flag step
    @param lambda_n:
    @param lambda_k:
    @param U_n:
    @param flag:
    '''
    U_k = U_n
    f_bar = 1.

    p_x = [U_n]
    p_y = [lambda_n]
    K = get_dif_K( lambda_n, U_n )
    h, w = get_dif_d( lambda_n, U_n )
    g = get_dissipation( lambda_n, lambda_k, U_n, U_k, flag )
    res = lambda_k * f_bar - fn.subs( x, U_k )

    if h == 0.:
        d_U = ( lambda_k - lambda_n ) * f_bar / K
        d_lambda = 0.
        a = 1.

    else:
        A = array( [[K, -f_bar], [h, w]], dtype = float )
        rhs = array( [res, -g + tau] )
        a = linalg.solve( A, rhs )
        d_U = a[0]
        d_lambda = a[1]
    U_k += d_U
    lambda_k += d_lambda
    p_x.append( U_k )
    p_y.append( lambda_k )

    k = 1
    while linalg.norm( a ) > tol:
        if flag == 'iter':
            K = get_dif_K( lambda_k, U_k )
            h, w = get_dif_d( lambda_k, U_k )
            A = array( [[K, -f_bar], [h, w]], dtype = float )
            #print "A ", A
        elif flag == 'K-iter':
            K = get_dif_K( lambda_k, U_k )
            A = array( [[K, -f_bar], [h, w]], dtype = float )
            #print "A ", A
        g = get_dissipation( lambda_n, lambda_k, U_n, U_k, flag )
        res = lambda_k * f_bar - fn.subs( x, U_k )
        rhs = array( [res, -g + tau] )
        #print "rhs ", rhs
        a = linalg.solve( A, rhs )
        d_U = a[0]
        U_k += d_U
        #print 'U_k ', U_k
        d_lambda = a[1]
        lambda_k += d_lambda
        #print 'lambda_k ', lambda_k
        p_x.append( U_k )
        p_y.append( lambda_k )
        k += 1
    print "n iter ", k
    #print "eq. state ", p_x[-1], " ", p_y[-1]
    pylab.plot( p_x, p_y , 'o-' )
    return p_x[-1], p_y[-1]

def get_val_numPOST( lambda_n, lambda_k, U_n, flag = 'step' ):
    '''
    dissipation control - d_lambda computed from d_U [Eckhard]
    @param lambda_n:
    @param lambda_k:
    @param U_n:
    @param flag:
    '''
    U_k = U_n
    f_bar = 1.

    p_x = [U_n]
    p_y = [lambda_n]
    K = get_dif_K( lambda_n, U_n )
    d_U = ( lambda_k - lambda_n ) * f_bar / K
    #U_k += d_U
    print 'd_U ', d_U
    #print 'U_k ', U_k

    while d_U > tol:
        d_lambda = ( 2 * tau - ( lambda_n * ( U_k + d_U - U_n ) - \
                                U_n * ( lambda_k - lambda_n ) ) * f_bar ) / \
                                 ( lambda_n * d_U - U_n )

        d_U = ( lambda_k + d_lambda - lambda_n ) * f_bar / K
        print 'd_U ', d_U
    pylab.plot( p_x, p_y, 'o-' )



if __name__ == '__main__':
    lambda_n = .5
    U_n = get_val_ana( lambda_n )
    print 'U_n ', U_n
    print '==============================='
    #compute the equlibrium with following methods
    #NR from U_n
    get_val_numNR( lambda_n, 1., U_n, 'iter' )
    #modified NR from U_n
    get_val_numNR( lambda_n, 1., U_n, 'step' )
    #dissipation control from U_n - linearized in step
    U_k, lambda_k = get_val_numD( lambda_n, lambda_n, U_n, 'step' )
    #TODO:dissipation control from origin - linearized in iteration
    #get_val_numD( 0.0, 0., 0.0, 'iter' )
    #dissipation control from U_n - linearized in iteration
    U_k, lambda_k = get_val_numD( lambda_n, lambda_n, U_n, 'iter' )
    #dissipation control from U_n - K linearized in iteration
    #U_k, lambda_k = get_val_numD( lambda_n, lambda_n, U_n, 'K-iter' )

    #compute the dissipation analytically
#    print 'numerical - analytical misfit %d'\
#            % ( ( compute_dissipation( lambda_n, lambda_k, U_n , U_k ) / \
#                  tau - 1. ) * 100. ), '%'
    #print compute_dissipation( 0., 1., 0. , 1. )

    plot_fn()
    pylab.ylim( 0., 1.1 )
    pylab.show()
