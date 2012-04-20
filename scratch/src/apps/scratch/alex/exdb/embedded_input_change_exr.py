'''
exTreated on exTMpr 30, 2010

@author: alexander
'''

from enthought.traits.api import \
    HasTraits, Instance, Int, on_trait_change, \
    Event, Enum, DelegatesTo, Property

import os

import pickle

from os.path import \
    exists




class CM( HasTraits ):
    
    value = Int( 0, input = True )

    input_change = Event
    @on_trait_change('+input')
    def _set_input_change(self):
        print 'CM: input change raised'
        self.input_change = True
    
    
class CCS( HasTraits ):

    rho_cc = Int(3, input = True)
    
    cm_key = Enum( 'key_a','key_b', input = True )

    cm = Instance( CM )
    def _cm_default(self):
        return CM()
    
    input_change = Event
    @on_trait_change('+input, cm.input_change')
    def _set_input_change(self):
        print 'CCS: input change propagated'
        self.input_change = True

class ExT( HasTraits ):
    
    
    ccs = Instance( CCS, input = True )
    def _ccs_default(self):
        return CCS()

    input_change = Event
    @on_trait_change('+input, ccs.input_change')
    def _set_input_change(self):
        print 'ExT: input change propagated'
        self.input_change = True

    rho_cc = DelegatesTo('ccs')
    

class ExR( HasTraits ):
    
    ext = Instance( ExT )
    def _ext_default(self):
        return ExT()
    
    input_change = Event
    @on_trait_change('+input, ext.input_change')
    def _set_input_change(self):
        print 'ExR: input change received'
        self.input_change = True


def delete_pickle_file( data_file ):
    print 'XXX Delete the pickle file if it exists--------------------'
    
    dir_path  = os.path.dirname( data_file )
    file_name = os.path.basename( data_file )
    file_split = file_name.split('.')
    pickle_file_name = os.path.join( dir_path, file_split[0] + '.pickle' ) 
    if os.path.exists( pickle_file_name ):
        print '--- pickle file removed: ', pickle_file_name, ' ---'
        os.remove( pickle_file_name )
    else:
        print '--- pickle file does not exist: ', pickle_file_name, ' ---'





if os.path.exists( 'ext.pickle' ):
    print 'pickle file removed'
    os.remove( 'ext.pickle' )
else:
    print 'pickle file does not exist'
    
print 'construct exr'
exr = ExR()        

print 'initial key'
print exr.ext.ccs.cm_key

print 'saving ext'
file_name = 'ext.pickle'
file = open( file_name,'w')
pickle.dump( exr.ext, file )
file.close()

print 'changing key'
exr.ext.ccs.cm_key = 'key_b'
print 'changed key'
print exr.ext.ccs.cm_key


