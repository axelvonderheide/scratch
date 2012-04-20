'''
Created on May 27, 2009

@author: jakub
'''
from enthought.traits.api import \
    HasTraits
    
class MyObj(HasTraits):
    
    def __getitem__(self, idx):
        
        print isinstance(idx,tuple)
        
if __name__ == '__main__':
    obj = MyObj()
    obj[0:1,2]