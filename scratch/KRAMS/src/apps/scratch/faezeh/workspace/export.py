from enthought.traits.api import Float, HasTraits, Button
from enthought.traits.ui.menu import OKButton, CancelButton
from enthought.traits.ui.api import View, Item, InstanceEditor
import os

class InputParam(HasTraits):

    height = Float()
    width  = Float()
    
    export_param = Button()
    def _export_param_fired(self):
        test_file = os.path.join('', 'Export_file')
        output_file = open(test_file + '.rtf','w')
        output_file.write(self.height.__str__() + '\n'+ self.width.__str__())
      
    view_traits = View( Item("height"),
                        Item("width"),
                        Item('export_param', label = 'Export to file', style='simple', show_label = False),
                        resizable = True,
                        buttons = [ OKButton, CancelButton ],
                        height = 0.5,
                        width = 0.5 )
    
if __name__ == '__main__': 
    ip = InputParam()
    ip.configure_traits()    
