'''
Created on May 29, 2009

@author: rchx
'''
from enthought.traits.api import HasTraits, Str, List, Float

class _LItem( HasTraits ):
    
    name = Str('')
    
    def export_latex( self, object ):
        print 'name', self.name, 'value', getattr( object, self.name )
        

def LItem( name, **args ):
    return _LItem( name = name, **args )

class LView( HasTraits ):
    '''
    classdocs
    '''
    
    # list of latex items
    #
    litems = List( _LItem )

    def __init__(self, *args, **kargs ):
        '''
        Constructor
        '''
        super( LView, self ).__init__(**kargs)
        for value in args:
            self.litems.append( value )
        
class MaterialModel( HasTraits ):
    
    E_modulus = Float( 30e3 )
    nu = Float( 0.2 )
    
    lview_trait = LView( LItem('E_modulus' ),
                         LItem('nu') )
    
mm = MaterialModel()

class LatexEngine( HasTraits ):
    
    def export_obj( self, obj, lview = 'lview_trait' ):
        
        if hasattr( obj, lview ):
            for litem in getattr( obj, lview ).litems:
                # now process the individual views
                litem.export_latex( obj )
        else:
            raise AttributeErrror, 'No latex view defined on %s' % obj


le = LatexEngine()
le.export_obj(mm)
    
    