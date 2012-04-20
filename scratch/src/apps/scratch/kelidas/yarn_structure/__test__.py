'''
Created on Dec 12, 2010

@author: kelidas
'''

import unittest
from numpy import min, array, all
from os.path import join
from promod.simdb import SimDB

from yarn_data_view import YarnRawData

simdb = SimDB()
test_dir = join( simdb.exdata_dir, 'trc', 'yarn_structure', 'SPEC' )
PX2MUM = 1.2

class TestYarnRawData( unittest.TestCase ):
    def setUp( self ):
        print 'setting up'
        self.yarn_raw_data = YarnRawData()
        self.yarn_raw_data.raw_data_directory = test_dir

    def test_cut_raw_data( self ):
        x_true = array( [[0, 0, 0, 0, 0, 0],
                         [0, 100, 50, 50, 200, 100],
                         [100, 300, 300, 400, 400, 350]] ).T / PX2MUM
        cf = array( [[0.25, 0., 0., 0., 0.55555556, 0.55555556],
                     [0.25, 0., 0.05555556, 0., 0.55555556, 0.55555556]] )

        cut_cfl_0 = array( [[0., 422.5, 422.5, 0., 0.],
                            [0., 0., 0., 0., 0.]] )
        real_cfl_0 = array( [[0., 489.91113955, 489.91113955, 0., 0.],
                             [0., 0., 0., 0., 0.]] )
        slack_0 = array( [[0., 0.13759871, 0.13759871, 0., 0],
                          [0., 0., 0., 0., 0.]] )
        
        cut_cfl_20 = array( [[0., 422.5, 422.5, 0., 0. ]
                             [0., 422.5, 422.5, 0., 0.]] )
        real_cfl_20 = array( [[0., 489.91113955, 489.91113955, 0., 0.]
                              [0., 527.78071613, 527.78071613, 0., 0.]] )
        slack_20 = array( [[0., 0.13759871, 0.13759871, 0., 0.]
                           [0., 0.19947814, 0.19947814, 0., 0.]] )
        
        self.assertEqual( all( self.yarn_raw_data.cut_raw_data[0] == x_true.flatten() ), True )
        self.assertEqual( all( self.yarn_raw_data.cut_raw_data[0] == cf.flatten() ), True )
        self.assertEqual( all( self.yarn_raw_data.cut_raw_data[0] == cut_cfl_0.flatten() ), True )
        self.assertEqual( all( self.yarn_raw_data.cut_raw_data[0] == real_cfl_0.flatten() ), True )
        self.assertEqual( all( self.yarn_raw_data.cut_raw_data[0] == slack_0.flatten() ), True )
        self.assertEqual( all( self.yarn_raw_data.cut_raw_data[0] == cut_cfl_20.flatten() ), True )
        self.assertEqual( all( self.yarn_raw_data.cut_raw_data[0] == real_cfl_20.flatten() ), True )
        self.assertEqual( all( self.yarn_raw_data.cut_raw_data[0] == slack_20.flatten() ), True )
        
        print 'cut_raw_data tested'


if __name__ == '__main__':
    unittest.main()

