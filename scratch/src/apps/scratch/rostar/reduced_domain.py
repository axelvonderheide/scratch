'''
Created on Oct 7, 2011

@author: rostar
'''

from iid_order_stats import IIDOrderStats
from mathkit.mfn.mfn_line.mfn_line import MFnLineArray
from enthought.traits.api import Property, cached_property
from stats.misc.generated_distribution import GenDistr
from scipy.optimize import broyden2 as slv
import numpy as np

class RD( IIDOrderStats ):

    pdfk = Property( depends_on = 'k, n, x_arr, distr' )
    @cached_property
    def _get_pdfk( self ):
        return MFnLineArray( xdata = self.x_arr, ydata = self.kth_pdf() )

    gd = GenDistr()

    def opt( self, pdf ):
        PDFK = self.pdfk.get_values( self.x_arr )
        self.gd.pdf_values = pdf
        self.gd.x_values = self.x_arr
        sf = self.gd.sf( self.x_arr )
        N = self.n - self.k + 1.
        K = 1.
        val = PDFK - N * np.array( pdf ) * sf ** ( N - K )
        res = val + np.abs( ( np.trapz( pdf, self.x_arr ) - 1 ) )
        return res

    reduced_rv = Property( depends_on = 'k, n, x_arr, distr' )
    @cached_property
    def _get_reduced_rv( self ):
        return slv( self.opt, self.distr.pdf( self.x_arr ) )

if __name__ == '__main__':

    from scipy.stats import norm
    import pylab as p

    rd = RD( distr = norm, k = 2, n = 10, x_arr = np.linspace( -6, 3, 2000 ) )
    amp_distr = rd.reduced_rv
    p.plot( rd.x_arr, amp_distr, label = 'AMP' )
    print 'AMP', np.trapz( amp_distr, rd.x_arr )
    p.plot( rd.x_arr, rd.kth_pdf(), label = 'kth', color = 'black' )
    print 'kth', np.trapz( rd.kth_pdf(), rd.x_arr )
    p.plot( rd.x_arr, rd.distr.pdf( rd.x_arr ), label = 'IID', ls = '--' )
    print 'IID', np.trapz( rd.kth_pdf(), rd.x_arr )
    iid = RD( distr = GenDistr( x_values = rd.x_arr, pdf_values = np.array( amp_distr ) ), x_arr = rd.x_arr, n = 9, k = 1 )
    p.plot( iid.x_arr, iid.kth_pdf(), label = 'kthIID', color = 'red', lw = 2, ls = '--' )
    p.legend()
    p.show()
