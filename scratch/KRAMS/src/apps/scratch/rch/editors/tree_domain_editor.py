#  Copyright (c) 2007, Enthought, Inc.
#  License: BSD Style.

# tree_editor.py -- Example of a tree editor

#--[Imports]--------------------------------------------------------------------
from enthought.traits.api \
    import HasTraits, Str, Regex, List, Instance
from enthought.traits.ui.api \
    import TreeEditor, TreeNode, View, Item, VSplit, \
           HGroup, Handler, Group
from enthought.traits.ui.menu \
    import Menu, Action, Separator
from enthought.traits.ui.wx.tree_editor \
    import NewAction, CopyAction, CutAction, \
           PasteAction, DeleteAction, RenameAction

#--[Code]-----------------------------------------------------------------------

# DATA CLASSES

class FERefinementLevelGrid ( HasTraits ):
    name  = Str( '<unknown>' )
    title = Str
    phone = Regex( regex = r'\d\d\d-\d\d\d\d' )
    
    def default_title ( self ):
        self.title = 'Senior Engineer'
    
class FEGrid ( HasTraits ):
    name      = Str( '<unknown>' )
    children = List( FERefinementLevelGrid )

class FEDomainList ( HasTraits ):
    name        = Str( '<unknown>' )
    subdomains = List( FEGrid )
          
class FEDomainTree ( HasTraits ):
    name    = Str( '<unknown>' )
    domain_list = Instance( FEDomainList )

# INSTANCES

jason = FERefinementLevelGrid( 
     name  = 'Jason',
     title = 'Engineer', 
     phone = '536-1057' )
     
mike = FERefinementLevelGrid( 
     name  = 'Mike',
     title = 'Sr. Marketing Analyst', 
     phone = '536-1057' )
     
dave = FERefinementLevelGrid(
     name  = 'Dave',
     title = 'Sr. Engineer',
     phone = '536-1057' )
     
susan = FERefinementLevelGrid(
     name  = 'Susan',
     title = 'Engineer',
     phone = '536-1057' )
     
betty = FERefinementLevelGrid(
     name  = 'Betty',
     title = 'Marketing Analyst' )
        
owner = FEDomainTree(
    name    = 'wile',
    domain_list = FEDomainList( 
        name = 'Domain Tree',
        subdomains = [
            FEGrid( 
                name = 'Grid 0',
                children = [ mike, betty ]
            ),
            FEGrid(
                name = 'Grid 1',
                children = [ dave, susan, jason ] 
            )
        ]
    )
)

print owner.domain_list.subdomains

# Actions used by tree editor context menu

# Tree editor 
tree_editor = TreeEditor( hide_root = True,
    nodes = [
        TreeNode( node_for  = [ FEDomainList ],
                  auto_open = True,
                  children  = 'subdomains',
                  label     = '=FEDomainList' ),
        TreeNode( node_for  = [ FEGrid ],
                  auto_open = True,
                  children  = 'children',
                  label     = '=FEGrid',
                  ),
        TreeNode( node_for  = [ FERefinementLevelGrid ],
                  auto_open = True,
                  label     = '=FERefinementLevelGrid'
                  )
    ]
)
# The main view
view = View( 
           Group( 
               Item( 
                    name = 'domain_list',
                    id = 'domain_list',
                    editor = tree_editor, 
                    resizable = True,
                    show_label = False ), 
                orientation = 'vertical',
                show_labels = True,
                show_left = True, ),
            title = 'FEDomainList Structure',
            id = \
             'enthought.traits.ui.tests.tree_editor_test',
            dock = 'horizontal',
            drop_class = HasTraits,
            buttons = [ 'Undo', 'OK', 'Cancel' ],
            resizable = True,
            width = .3,
            height = .3 )
                       
if __name__ == '__main__':
    owner.configure_traits( view = view )

