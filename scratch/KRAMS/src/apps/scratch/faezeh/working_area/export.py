from enthought.traits.api import Float, HasTraits, Button
from enthought.traits.ui.menu import OKButton, CancelButton
from enthought.traits.ui.api import View, Item, InstanceEditor
import os

class Exportable( HasTraits ):
    
    eparams = []
    
    def export(self, file):
        eparams = self.eparams
        for e in eparams:
            file.write( e + ' = ' + str( getattr( self, e ) ) + '\n' )
            if isinstance( e, Exportable ):
                self.export( file )

class InputParam( Exportable ):

    height = Float()
    width  = Float()
    
    eparams = ["height", "width" ]
    
    view_traits = View( Item("height"),
                        Item("width"),
                        Item('export_param', label = 'Export to file', style='simple', show_label = False),
                        resizable = True,
                        buttons = [ OKButton, CancelButton ],
                        height = 0.5,
                        width = 0.5 )
    
def export_rtf( exported_obj ):
    
    # get the parameters to be exported
    
    test_file = os.path.join('', 'Export_file')
    output_file = open(test_file + '.rtf','w')

    exported_obj.export( output_file )

    output_file.close()

ip = InputParam()

export_rtf( ip )