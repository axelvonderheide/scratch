'''
Created on Apr 28, 2010

@author: alexander
'''
from enthought.traits.api import \
    HasTraits, Directory, List, Int, Float, Any, \
    on_trait_change, File, Constant, Instance, Trait, \
    Array, Str, Property, cached_property, WeakRef, \
    Dict, Button, Bool, Enum, Event, implements, \
    DelegatesTo, Regex

from enthought.traits.ui.api import \
    View, Item, DirectoryEditor, TabularEditor, HSplit, VGroup, \
    TableEditor, EnumEditor, Handler, FileEditor, VSplit, Group, \
    InstanceEditor, HGroup, Spring
    

from enthought.units.ui.quantity_view import \
    QuantityView
    
from enthought.units import \
    *
    
from enthought.units.convert import \
    convert 


from enthought.units.quantity import \
    Quantity    

from scipy import arange

# Major package imports
#from enthought.util.numerix import array, arange
#from enthought.util.testingx import assert_equal, assert_almost_equal

# Local units imports
import enthought.units as units
from enthought.units.mass import kg, metric_ton
from enthought.units.temperature import *
from enthought.units import area, density, energy, force, length, mass, power, \
              pressure, speed, substance, time, frequency, acceleration, \
              temperature
from enthought.units import unit_system
from enthought.units.quantity import Quantity
from enthought.units.style_manager import style_manager
from enthought.units.unit_manager import unit_manager
from enthought.units.smart_unit import is_dimensionless
    
from enthought.units.SI import dimensionless
from enthought.units.speed import meters_per_second
    
from enthought.units.unit_parser import unit_parser


## scalar example
from enthought.units.SI \
    import meter, second

from enthought.units.length \
    import *
    


## an array example
from scipy import arange


class MetaQuantityView( View ):
    """ Default Traits View for MetaQuantity objects. """
    
    def __init__(self, *args, **traits):
        """ Create a new MetaQuantityView. """
        
        handler = traits.setdefault('handler', MetaQuantityViewHandler())
        handler.known_names = traits.pop('known_names', [])
        handler.any_name    = traits.pop('any_name', True)
        
        if handler.any_name:
            evaluate = handler.validate_name
        else:
            evaluate = None
        
        name_editor = EnumEditor( name="known_names", object='handler',
                                  evaluate=evaluate)

        name_item = Item( name='name', label='Name',
                                editor=name_editor, id='name_item' )

        super( MetaQuantityView, self ).__init__(
            Item( name='name', editor=name_editor),
            Item( name='family_name', label='Measure of' ),
            Item( name='units' ),
            *args,
            **traits
            )


class MetaQuantityViewHandler(Handler):
    """ The MetaQuantityViewHandler manages the optional limiting of name
    to selection from a predefined list.
    """
    
    # User is limited to selecting the name of a known Quantities.
    known_names = List
    
    # When True, user is not limited to the known names list.
    any_name = Bool(True)
   
    def validate_name(self, name):
        """ Validate name against the known names and any_name flag. 
        Returns name if validation passes. 
        Raises TraitError otherwise.
        """
        name = name.strip()
        
        if self.any_name and len(name) == 0:
            raise TraitError, 'name must be specified'
        
        if not (self.any_name or name in self.known_names):
            raise TraitError, 'invalid name %s' % name
        
        return name

# ------------------------------------------------------------------

class A( HasTraits ):
    '''Example class
    '''
    known_names = ['Q','D']
    b = Int(3, unit = 'cm')
    c = Int(3 )
    Q = Quantity(2.2, meter/second)
        
        
    ## an array example

    data = arange(3.)
    D = Quantity(data, meter, 'depth')

#    traits_view = MetaQuantityView( 
#    traits_view = QuantityView(
    traits_view = View(
##                       Item('b', show_label = True,),
##                       Item('c', show_label = True,),
                       Item('Q', show_label = True,),
                       Item('D', show_label = True,),
                       width = 0.4,
                       height = 0.4,
                       scrollable = True,
                       resizable = True,
                       )

#---------------------------------

if __name__ == '__main__':

    Q = Quantity(2.2, meter/second)
    #QQ =
    print 'Q',Q*Q
    data = arange(3.)
    print Quantity(data, meter, 'depth')


    A = A()
    A.configure_traits()
