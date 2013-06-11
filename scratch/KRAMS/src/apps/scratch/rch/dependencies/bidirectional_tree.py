
from enthought.traits.api import HasTraits, List, Str, Instance, This, WeakRef, on_trait_change

from enthought.traits.ui.api import View, Item, Group, TreeEditor, TreeNode, Handler

from enthought.traits.ui.menu \
    import Menu, Action, Separator
from enthought.traits.ui.wx.tree_editor \
    import NewAction, CopyAction, CutAction, \
           PasteAction, DeleteAction, RenameAction

add_node_action = Action(
    name='Add',
    action='handler.add_node(editor,object)')

get_root_action = Action(
    name='Get root',
    action='object.get_root().configure_traits()'
    #action='handler.get_root(editor,object)'
    )

class TreeHandler ( Handler ):
    
    def add_node ( self, editor, object ):
        object.subnodes.append( Node( ) )            

    def get_root ( self, editor, object ):
        object.get_root().configure_traits( kind = 'modal' ) 
                   
class Node( HasTraits ):
    
    subnodes = List

    name = Str( auto_set = False, enter_set = True )

    parent = WeakRef
    
    def get_root(self):
        print 'parent'
        print self.parent
        if self.parent == None:
            return self
        else:
            return self.parent.get_root()
    
    @on_trait_change('subnodes.+')
    def _reset_parent(self):
        for n in self.subnodes:
            n.parent = self
    
    traits_view = View( Item('name'), 
                        Item( 'parent@' ),
                        Item('subnodes@'),
                        resizable = True,
                        scrollable = True,
                        width = 200, height = 200 )

tree_editor = TreeEditor(
    nodes = [
        TreeNode( node_for = [ Node ],
                  auto_open = True,
                  label = 'name',
                  children = 'subnodes',
                  menu      = Menu( NewAction,
                                    Separator(),
                                    DeleteAction,
                                    Separator(),
                                    RenameAction,
                                    Separator(),
                                    CopyAction, 
                                    CutAction, 
                                    PasteAction,
                                    Separator(),
                                    get_root_action ),    
                   add = [Node],
                   ),
        ],
    hide_root = False
    )

class Tree( HasTraits ):
    
    root = Instance( Node )
    
    traits_view = View( Item( name = 'root',
                              id = 'domain_list',
                              editor = tree_editor, 
                              show_label = False ),
                              id = 'simvisage.ibvpy.mesh.fe_domain_tree',
                              dock = 'horizontal',
                              handler = TreeHandler(),
                              resizable = True,
                              scrollable = True )

if __name__ == '__main__':
    
    r = Node( name = 'n1', subnodes = [Node( name = 'n2', subnodes = [ Node( name = 'n3' ) ] ),
                                       Node( name = 'n4', subnodes = [ Node( name = 'n5' ) ] ) ] )
    
    t = Tree( root = r )
    t.configure_traits()
    
