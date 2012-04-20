'''
Created on 31.8.2011

class for evaluating the kth statistics of a subset of n realizations
all of which have a common distribution.
Example: k = 1, n = 10 evaluates the PDF of the weakest of 10 elements

@author: Q
'''

from enthought.traits.api import HasTraits, Int

import scipy as sp
import numpy as np

class IIDOrderStats( HasTraits ):

    distr = 'distribution'

    n = Int( 10, auto_set = False, enter_set = True,
            desc = 'number of realizations' )

    k = Int( 1, auto_set = False, enter_set = True,
            desc = 'kth order statistics to evaluate' )

    x_arr = np.linspace( 0, 13, 300 )

    def n_realizations( self ):
        return self.distr.ppf( np.random.rand( self.n ) )

    def kth_pdf( self ):
        '''evaluates the PDF of the kth entry'''
        n = self.n
        k = self.k
        x = self.x_arr
        fct = sp.misc.factorial
        pdf = self.distr.pdf
        cdf = self.distr.cdf
        constant = fct( n ) / ( fct( k - 1 ) * fct( n - k ) )
        return constant * pdf( x ) * cdf( x ) ** ( k - 1 ) * ( 1 - cdf( x ) ) ** ( n - k )

if __name__ == '__main__':

    from matplotlib import pyplot as plt
    from scipy.stats import norm

    n2_list = []
    n7_list = []
    n12_list = []
    for i in range( 50000 ):
        n = norm( 6, 2 ).ppf( np.random.rand( 12 ) )
        n.sort()
        n2_list.append( n[1] )
        n7_list.append( n[6] )
        n12_list.append( n[11] )


    iid = IIDOrderStats( distr = norm( 6, 2 ), n = 12, k = 2 )
    pdf1 = iid.kth_pdf()
    iid.k = 7
    pdf2 = iid.kth_pdf()
    iid.k = 12
    pdf3 = iid.kth_pdf()
    plt.hist( n2_list, 50, normed = True, label = 'monte carlo k = 2' )
    plt.hist( n7_list, 50, normed = True, label = 'monte carlo k = 7' )
    plt.hist( n12_list, 50, normed = True, label = 'monte carlo k = 12' )
    plt.plot( iid.x_arr, pdf1, color = 'black', lw = 2 )
    plt.plot( iid.x_arr, pdf2, color = 'black', lw = 2 )
    plt.plot( iid.x_arr, pdf3, color = 'black', lw = 2, label = 'prediction k = 2,7,12' )
    plt.title( 'comparsion of the implemented formula with an MC simulation' )
    plt.legend()
    plt.show()
