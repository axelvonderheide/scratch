'''
Created on Apr 1, 2010

@author: alexander
'''
from enthought.traits.api import \
    HasTraits, Directory, List, Int, Float, Any, \
    on_trait_change, File, Constant, Instance, Trait, \
    Array, Str, Property, cached_property, WeakRef, \
    Dict, Button, Bool, Enum, Event, implements, \
    DelegatesTo, Callable
    
from numpy import \
    log
    
from promod.simdb.simdb_class import \
    SimDBClass, SimDBClassExt
    
class ConcreteMixture( SimDBClass ):
    '''Describes the properties of the concrete matrix
    '''
    
    # E-modulus of the concrete after 28d
    E_m28 = Float( unit = 'MPa', simdb = True, input = False, auto_set = False, enter_set = False )

    # developement of the E-modulus depending on the age at the time of testing:
    get_E_m_time = Callable
    
# Setup the database class extension 
#
ConcreteMixture.db = SimDBClassExt(
            klass = ConcreteMixture,
            constants = { 
                'PZ-0708-1' : ConcreteMixture(
                                           E_m28 = 33036,
                                           get_E_m_time = lambda t: 4665. * log( t + 0.024 ) + 17487. 
                                           ),
                'FIL-10-09' : ConcreteMixture(
                                           E_m28 = 28700,
                                           # function for the evolution derived based on only 
                                           # three values: Em0 = 0, Em7 = 23600, Em28 = 28700
                                           get_E_m_time = lambda t: 3682. * log( t + 0.012 ) + 16429.  
                                           )}
            )

if __name__ == '__main__':
    ConcreteMixture.db.configure_traits()