
    
from enthought.traits.ui.wx.editor \
    import Editor
    
import wx

class _FooEditor( Editor ):
    def init(self, parent):
        factory = self.factory
        self.control = self._create_canvas( parent )

    def _create_canvas( self, parent ):
        panel = wx.Panel(parent, -1, style=wx.CLIP_CHILDREN)
        sizer = wx.BoxSizer(wx.VERTICAL)
        panel.SetSizer(sizer)

        print 'a', self.context_object.a
        print 'b', self.context_object.b
        
        print 'name', self.name
        print 'value', self.value

        return panel
        
    def update_editor(self):
        print 'Im updating'
        pass

from enthought.traits.ui.basic_editor_factory \
    import BasicEditorFactory

class FooEditor( BasicEditorFactory ):
    
    klass = _FooEditor