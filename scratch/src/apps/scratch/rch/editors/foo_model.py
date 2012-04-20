
from enthought.traits.api import HasTraits, Float, Instance
from enthought.traits.ui.api  import View, Item
from foo_editor import FooEditor

class FooModel(HasTraits):
    
    a = Float(10)
    b = Float(20)
    
if __name__ == '__main__':
    
    class FooContainer( HasTraits ):
        foo_trait = Instance( FooModel )
    
    fm = FooModel()
    fm.configure_traits( view = View( Item('foo_trait@s', editor = FooEditor() ) ) )