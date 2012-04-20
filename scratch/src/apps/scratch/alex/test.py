#### Deprecated file - to be removed [rch]

from enthought.traits.api import \
     Array, Bool, Enum, Float, HasTraits, HasStrictTraits, \
     Instance, Int, Trait, Str, Enum, \
     Callable, List, TraitDict, Any, Range, \
     Delegate, Event, on_trait_change, Button, \
     Interface, WeakRef, implements, Property, cached_property, Tuple, \
     Dict
from enthought.traits.ui.api import Item, View, HGroup, ListEditor, VGroup, \
     HSplit, Group, Handler, VSplit, TableEditor, ListEditor

from enthought.traits.ui.menu import NoButtons, OKButton, CancelButton, \
     Action

from enthought.traits.ui.ui_editors.array_view_editor \
    import ArrayViewEditor

from enthought.traits.ui.table_column \
    import ObjectColumn, ExpressionColumn

from enthought.traits.ui.table_filter \
    import TableFilter, RuleTableFilter, RuleFilterTemplate, \
           MenuFilterTemplate, EvalFilterTemplate, EvalTableFilter

from numpy import ix_, mgrid, array, arange, c_, newaxis, setdiff1d

from core.rv import RTrace

# tvtk related imports
#
from enthought.traits.ui.api import View, Item, HSplit, VSplit, InstanceEditor
from enthought.tvtk.api import tvtk
from core.sctx import \
    ISDomain, SDomain

from enthought.pyface.tvtk.scene_editor import SceneEditor




class Parent(HasTraits):
    last_name = Str('')
    

class Child(HasTraits):
    age = Int
    father = Instance(Parent)
    last_name = Delegate('father')
    def _age_changed(self,old,new):
        print 'Age changed from %s to %s' %(old, new)
        
joe = Parent()
joe.last_name = 'Johnson'
moe = Child()
moe.father = joe

moe.configure_traits()




