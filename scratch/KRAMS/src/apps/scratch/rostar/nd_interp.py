'''
Created on Sep 22, 2011

@author: rostar
'''

from enthought.traits.api import HasTraits, Array, List
import numpy as np
from scipy import ndimage

class NDIdxInterp( HasTraits ):

    data = Array # nd array of values (measured, computed..) of size orthogonalize(axes_values)
    axes_values = List( Array( float ) ) # list of control input parameter values

    def __call__( self, *gcoords, **kw ):
        ''' 
        kw: dictionary of values to interpolate for;
        len(kw) has to be equal the data dimension
        '''
        order = kw.get( 'order', 1 )
        if len( self.axes_values ) != len( gcoords ):
            raise TypeError, 'method takes {req} arguments ({given} given)'.format( req = \
                len( self.axes_values ), given = len( gcoords ) )

        # the value to be interpolated is defined by coordinates assigning
        # the data grid; e.g. if coords = array([[1.5, 2.5]]) then coords.T
        # point to the half distance between the 1. and 2. entry in the first
        # dimension and the half distance between the 2. and 3. entry of the 2. dim
        icoords = self.get_icoords( gcoords )
        # specify the data for the interpolation
        data = self.data
        # interpolate the value (linear)
#        a, b, c = [0.5, 0.5, 0.5], [0, 0, 0], [0, 1, 2]
#        icoords = [a, b, c]
        val = ndimage.map_coordinates( data, icoords, order = order, mode = 'nearest' )
        return val.flatten()

    def get_icoords( self, gcoords ):
        '''
        gcoords: values to be interpolated for
        this method transforms the global coords to "index" coords
        '''
        icoords = [np.interp( gcoord, axis_values, np.arange( len( axis_values ) ) )
                   for gcoord, axis_values in zip( gcoords, self.axes_values ) ]
        # ndimage.map_coordinates takes the reverse order of the list
        icoords.reverse()
        return icoords

if __name__ == '__main__':

    data = np.arange( 27 ).reshape( 3, 3, 3 )
    print 'data =', data

    axes_values = [[0.1, 0.2, 0.3], [10, 20, 30], [0, 1, 2]]
    ni = NDIdxInterp( data = data, axes_values = axes_values )
    w, x, f = np.linspace( 0.1, 0.3, 3 ), np.linspace( 10, 30, 3 ), np.linspace( 0, 2, 3 )

    res = ni( w, x, f )
    print 'result', res
