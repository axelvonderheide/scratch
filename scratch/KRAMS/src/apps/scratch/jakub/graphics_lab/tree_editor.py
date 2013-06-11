'''
Created on May 6, 2009

@author: jakub
'''
#  Copyright (c) 2007, Enthought, Inc.
#  License: BSD Style.#-- Imports --------------------------------------------------------------------

from os \
    import getcwd
    
from enthought.traits.api \
    import HasTraits, Property, Directory, adapts, property_depends_on
    
from enthought.traits.ui.api \
    import View, VGroup, Item, TreeEditor, ITreeNode, ITreeNodeAdapter
    
from enthought.io.api \
    import File
    
#-- FileAdapter Class ----------------------------------------------------------

class FileAdapter ( ITreeNodeAdapter ):
    
    adapts( File, ITreeNode )
    
    #-- ITreeNodeAdapter Method Overrides --------------------------------------

    def allows_children ( self ):
        """ Returns whether this object can have children.
        """
        return self.adaptee.is_folder

    def has_children ( self ):
        """ Returns whether the object has children.
        """
        children = self.adaptee.children
        return ((children is not None) and (len( children ) > 0))

    def get_children ( self ):
        """ Gets the object's children.
        """
        return self.adaptee.children
        
    def get_label ( self ):
        """ Gets the label to display for a specified object.
        """
        return self.adaptee.name + self.adaptee.ext
        
    def get_tooltip ( self ):
        """ Gets the tooltip to display for a specified object.
        """
        return self.adaptee.absolute_path
        
    def get_icon ( self, is_expanded ):
        """ Returns the icon for a specified object.
        """
        if self.adaptee.is_file:
            return '<item>'
            
        if is_expanded:
            return '<open>'
            
        return '<open>'

    def can_auto_close ( self ):
        """ Returns whether the object's children should be automatically 
            closed.
        """
        return True

#-- FileTreeDemo Class ---------------------------------------------------------

class FileTreeDemo ( HasTraits ):
    
    # The path to the file tree root:
    root_path = Directory( entries = 10 ) 
    
    # The root of the file tree:
    root = Property
    
    # The traits view to display:
    view = View(
        VGroup(
            Item( 'root_path' ),
            Item( 'root', 
                  editor = TreeEditor( editable = False, auto_open = 1 )
            ),
            show_labels = False
        ),
        width     = 0.33,
        height    = 0.50,
        resizable = True
    )
    
    #-- Traits Default Value Methods -------------------------------------------
    
    def _root_path_default ( self ):
        return getcwd()
    
    #-- Property Implementations -----------------------------------------------
    
    @property_depends_on( 'root_path' )
    def _get_root ( self ):
        return File( path = self.root_path )
    
#-- Create and run the demo ----------------------------------------------------

demo = FileTreeDemo()

# Run the demo (if invoked form the command line):
if __name__ == '__main__':
    demo.configure_traits()