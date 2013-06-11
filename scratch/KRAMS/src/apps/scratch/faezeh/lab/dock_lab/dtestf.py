from enthought.traits.api \
import *
#import HasTraits, Int, Str, View, Item
from enthought.traits.ui.api \
    import *
    
from enthought.traits.ui.menu \
    import *
           
from enthought.pyface.image_resource \
    import ImageResource

# Define the InputNumber class:
class InputNumber ( HasTraits ):
# Define the input value:
    value = Int( connect = 'from:the integer value' )
# Define the application object view:
    traits_view = View( 'value' )
# Define the ConvertNumber class:
class ConvertNumber ( HasTraits ):
# Define the number and converted value traits:
    number = Int( connect = 'to:the integer value' )
    text = Str
# Define the application object view:
    traits_view = View( Item( 'text', style='readonly' ) )
# Handle the 'number' trait being changed:
    def _number_changed ( self, number ):
        self.text = self._convert( number )
# Convert an integer to its text representation:
    def _convert ( self, n ):
        if n == 0:
            return 'zero'
            result = ''
        if n < 0:
            result = 'minus '
            n = -n

            if n >= 1000000000:
                result += \
                self._convert( n/1000000000 ) + ' billion '
                n %= 1000000000
            if n >= 1000000:
                result += \
                self._convert( n / 1000000 ) + ' million '
                n %= 1000000
            if n >= 1000:
                result += \
                self._convert( n / 1000 ) + ' thousand '
                n %= 1000
            if n >= 100:
                result += \
                self._convert( n / 100 ) + ' hundred '
                n %= 100
            if n >= 20:
                result += tens[ (n / 10) - 2 ]
                n %= 10
            if n > 0:
                result += '-'
            if n > 0:
                result += digits[ n ]
                return result.strip()
            
# List of multiples of ten names:
tens = [ 'twenty', 'thirty', 'forty', 'fifty',
'sixty', 'seventy', 'eighty', 'ninety' ]
# List of digit names:
digits = [ 'zero', 'one', 'two', 'three', 'four', 'five',
'six', 'seven', 'eight', 'nine', 'ten',
'eleven', 'twelve', 'thirteen', 'fourteen',
'fifteeen', 'sixteen', 'seventeen', 'eighteen',
'nineteen' ]
# Create importable application objects:
input_number = InputNumber()
convert_number = ConvertNumber()

if __name__ == '__main__':
    convert_number.configure_traits()
