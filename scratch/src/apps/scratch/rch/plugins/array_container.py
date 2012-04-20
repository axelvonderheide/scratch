

from enthought.traits.api import \
    Float, HasTraits, Int, List, Array, Dict, String, WeakRef, Instance, Enum, Property, Any, Delegate, \
    cached_property, Constant, on_trait_change

from enthought.traits.ui.api \
    import TabularEditor, View, Item, Group, HGroup, VGroup, HSplit, VSplit 
from enthought.traits.ui.menu \
    import OKButton, CancelButton
from enthought.traits.ui.tabular_adapter \
    import TabularAdapter


class ViewElem(HasTraits):
    '''Get the element numbers.
    '''
    elem_num = Int(-1)
    view = View( Item( 'elem_num', style = 'readonly' ) )


#----------------------------------------------------------------------------------
#                                     RIDVariable
#----------------------------------------------------------------------------------
class RIDVariable(HasTraits):
    """
    Association between a random variable and distribution.
    """
    varname = String

    traits_view = View( Group( Item( 'varname', style = 'readonly', show_label = False ) ),
                        resizable = True,
                        height=800 )

#-- Tabular Adapter Definition -------------------------------------------------

class RVAdapter ( TabularAdapter ):

    columns = [ ( 'Name',         'varname' ) ]
                
    font                      = 'Courier 10'
    variable_alignment        = Constant( 'right' )
    
#-- Tabular Editor Definition --------------------------------------------------

rtrace_editor = TabularEditor(
    selected   = 'current_variable',
    adapter    = RVAdapter(),
    operations = [ 'move' ],
#    auto_update = True
)

from array_object import ShowArray, tabular_editor
from numpy.random import random
#----------------------------------------------------------------------------------
#                                     SPIRRID                                     
#----------------------------------------------------------------------------------
class SPIRRID( HasTraits ):

    data = Array( Int, value = random( ( 1000, 3 ) ) )

    current_row = Int(-1)
    
    view_elem = Instance( ViewElem )
    def _view_elem_default(self):
        return ViewElem()
     
    @on_trait_change('current_row')
    def _update_view_elem(self):
        if self.current_row != -1:
            self.view_elem.elem_num = self.current_row
    
    @on_trait_change('current_row')
    def _display_current_row(self):
        print 'row value'
        print self.current_row
    
    variables = List
    def _variables_default( self ):
        '''
        reset the RIDVariable list
        '''
        varset = []
        for v in ['one','two','three']:
            varset.append (
                    RIDVariable( varname = v ) )
        return varset

    current_variable = Instance( RIDVariable )
    def _current_variable_default(self):
        return self.variables[0]
 
    @on_trait_change('current_variable')
    def _display_current_variable(self):
        print 'row value'
        print self.current_variable
           
    traits_view = View( VGroup( Item('data',editor = tabular_editor, style='readonly'),
                                Item('current_row'),
                                Item('view_elem@'),
                                VSplit( 
                                   Item('variables', show_label = False, editor = rtrace_editor ),
                                   Item('current_variable', show_label = False, style = 'custom', resizable = True),
                                   ),
                               ),
                        resizable = True,
                        width = 0.8,
                        height = 0.8,
                        )

if __name__ == '__main__':
    spirrid = SPIRRID()
    spirrid.configure_traits()