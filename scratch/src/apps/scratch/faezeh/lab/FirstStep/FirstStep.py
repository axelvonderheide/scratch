
# Imports:
from enthought.traits.api \
    import HasTraits
    
from math import cos

class X2COSX(HasTraits):
    ''' Simple class emulating a function x**2 * cos(x)
    '''
    
    def __call__(self, x):
        return x**2 * cos(x)
    
if __name__ == '__main__':
    fx = X2COSX()
    v = fx(2.1)
    print 'result', v
