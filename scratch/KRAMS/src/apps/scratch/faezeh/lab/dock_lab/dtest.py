from enthought.traits.api \
    import *
from enthought.traits.ui.api import *    
from enthought.traits.ui.table_column \
    import *
from enthought.developer.api \
    import *
# Define a read-only table column:
class ReadOnlyColumn ( ObjectColumn ):
    editable = False
# The table editor for displaying the FilePosition columns:
    files_table_editor = TableEditor(
    columns = [ ReadOnlyColumn( name = 'file_name' ),
                ReadOnlyColumn( name = 'name' ),
                ReadOnlyColumn( name = 'line' ),
                ReadOnlyColumn( name = 'column' ),
                ReadOnlyColumn( name = 'lines' ) ]
               )
# Define the application object class:
class FileDropper ( HasTraits ):
# The list of files being displayed:
    files = List( FilePosition, drop_file = True )
    # The table view for displaying them:
    traits_view = View(
                       Item( 'files',
                             show_label = False,
                             editor = files_table_editor
                             )
                       )
# Create an importable application object:
file_dropper = FileDropper()

if __name__ == '__main__':
    file_dropper.configure_traits()
