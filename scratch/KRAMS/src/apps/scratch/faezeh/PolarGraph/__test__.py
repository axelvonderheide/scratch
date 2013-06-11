
from mfn_polar import MFnPolar
import unittest

class TestSequenceFunctions(unittest.TestCase):
    
    def setUp(self):
        self.mp = MFnPolar( alpha = 0.5 )

    def test_get_value(self):
        '''
        make sure that the values for corner nodes get returned properly
        - testing (x,y) plane
        '''
        value = self.mp( 0.54 )
        self.assertEqual( value, 0.48 )        
