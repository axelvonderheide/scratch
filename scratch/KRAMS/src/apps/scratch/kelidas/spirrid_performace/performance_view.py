'''
Created on Jan 21, 2011

@author: kelidas
'''
from enthought.traits.api import HasTraits, Property, cached_property, \
    Event, Instance, File, Int, Range, \
     on_trait_change, Bool, Trait, Constant, Str
from enthought.traits.ui.api import Item, View, Group, \
    HGroup
from matplotlib.figure import Figure
from numpy import loadtxt, sum, max, min
from traits.editors.mpl_figure_editor import MPLFigureEditor



class Data( HasTraits ):
    '''
        data1 -- influence on compilation time
        data2 -- without influence on compilation time
    '''
    file1 = File( 'results/pullout_comb_int30_eps03_1.dat', filter = ['*.dat'] )
    file2 = Property( Str, depends_on = 'file1' )
    def _get_file2( self ):
        return self.file1[:-5] + '2.dat'

    data1 = Property( depends_on = 'file1' )
    def _get_data1( self ):
        data = loadtxt( self.file1, usecols = ( 1, 2, 5, 3, 4 ) )
        return data

    data1_text = Property( modified = True )
    def _get_data1_text( self ):
        return loadtxt( self.file1, usecols = [1], dtype = '|S1000', delimiter = '(' )

    data2 = Property( depends_on = 'file2' )
    def _get_data2( self ):
        data = loadtxt( self.file2, usecols = ( 1, 2, 5, 3, 4 ) )
        return  data

    data2_text = Property( modified = True )
    def _get_data2_text( self ):
        return loadtxt( self.file2, usecols = [1], dtype = '|S100', delimiter = '(' )

    traits_view = View( Item( 'file1', show_label = False ),
                        Item( 'file2', style = 'readonly', show_label = False ) )

class PView( HasTraits ):
    data_class = Instance( Data, () )

    data_sel = Trait( 'data2 -- no compilation time',
                     {'data1 -- compilation time':'data1',
                       'data2 -- no compilation time':'data2'}, modified = True )

    data = Property( depends_on = 'data_class' )
    def _get_data( self ):
        return getattr( self.data_class, self.data_sel_ )

    data_text = Property( depends_on = 'data_class' )
    def _get_data_text( self ):
        return getattr( self.data_class, self.data_sel_ + '_text' )

    zero = Constant( 0 )
    min_comb_slide = Property( Int, depends_on = 'comb_slider', modified = True )
    @cached_property
    def _get_min_comb_slide( self ):
        return int( min( self.data[:, 0] ) )
    max_comb_slide = Property( Int, modified = True, depends_on = 'comb_slider' )
    @cached_property
    def _get_max_comb_slide( self ):
        return int( max( self.data[:, 0] ) )

    comb_slider = Range( 'min_comb_slide', 'max_comb_slide', modified = True )

    max_slide = Property( Int, depends_on = 'comb_slider', modified = True )
    @cached_property
    def _get_max_slide( self ):
        return int( sum( self.data[:, 0] == self.comb_slider ) - 1 )
    sub_slider = Range( 'zero', 'max_slide', mode = 'slider', modified = True )
    sub_slider_on = Bool( False, modified = True )


    figure = Instance( Figure )
    def _figure_default( self ):
        figure = Figure()
        figure.add_axes( [0.1, 0.1, 0.8, 0.8] )
        return figure

    def plot_pullout( self ):
        ''' Plot pullout_comb_int30_eps03_1.dat
        '''
        import matplotlib.pyplot as plt
        from matplotlib import rc
        rc( 'font', family = 'serif',
             style = 'normal', variant = 'normal', stretch = 'normal', size = 16 )

        data = self.data_class.data2.copy()
        data[data < 0] = self.data_class.data1[data < 0]
        x = [1, 2, 3, 4]
        plt.figure( 0 )
        plt.subplots_adjust( wspace = 0.0, hspace = 0.0, bottom = .11 )#.5
        plt.plot( x, data[:, 1:][ data[:, 0] == 5].T, '-x', color = 'gray' )
        plt.xlim( 0.5, 4.5 )
        plt.ylim( 0, 1.1 * max( data[:, 1:][ data[:, 0] == 5] ) )

        plt.xlabel( '$\mathrm{configuration}$', size = 20 )
        plt.xticks( ( 1, 2, 3, 4 ) ,
                     ( '$\mathrm{I}$', '$\mathrm{II}$', '$\mathrm{III}$',
                       '$\mathrm{IV}$' ), size = 20, position = ( 0, -.01 ) )
        plt.ylabel( '$\mathrm{time\, [sec]}$', size = 20 )
        plt.yticks( ( 0, 100, 200, 300, 400, 500 ) ,
                     ( '$0$', '$100$', '$200$', '$300$', '$400$', '$500$' ) )
        plt.ylim( 0, 500 )
        tsize = 20
        #plt.text( 3.4, 60, 'D', size = tsize,
        #           horizontalalignment = 'center', verticalalignment = 'center' )
        #plt.text( 3.4, 150, 'C', size = tsize,
        #           horizontalalignment = 'center', verticalalignment = 'center' )
        #plt.text( 3.4, 300, 'B', size = tsize,
        #           horizontalalignment = 'center', verticalalignment = 'center' )
        #plt.text( 3.4, 400, 'A', size = tsize,
        #           horizontalalignment = 'center', verticalalignment = 'center' )
        plt.annotate( '$\mathcal{D}$', xy = ( 3.5, 60 ), xycoords = 'data',
                xytext = ( 3.8, 15 ), #textcoords = 'offset points',
                arrowprops = dict( arrowstyle = "->" ), size = tsize,
                )
        plt.annotate( '$\mathcal{C}$', xy = ( 3.4, 150 ), xycoords = 'data',
                xytext = ( 3.1, 150 ), #textcoords = 'offset points',
                arrowprops = dict( arrowstyle = "->" ), size = tsize,
                )
        plt.annotate( '$\mathcal{B}$', xy = ( 3.4, 270 ), xycoords = 'data',
                xytext = ( 3.65, 270 ), #textcoords = 'offset points',
                arrowprops = dict( arrowstyle = "->" ), size = tsize,
                )
        plt.annotate( '$\mathcal{A}$', xy = ( 3.5, 370 ), xycoords = 'data',
                xytext = ( 3.8, 410 ), #textcoords = 'offset points',
                arrowprops = dict( arrowstyle = "->" ), size = tsize,
                )
        plt.show()

    def plot_pullout_var2( self ):
        ''' Plot pullout_comb_int2000_eps03_var2_1.dat
        '''
        import matplotlib.pyplot as plt
        from matplotlib import rc
        rc( 'font', family = 'serif',
             style = 'normal', variant = 'normal', stretch = 'normal', size = 16 )

        data = self.data_class.data2.copy()
        data[data < 0] = self.data_class.data1[data < 0]
        x = [1, 2, 3, 4]
        plt.figure( 0 )
        plt.subplots_adjust( wspace = 0.0, hspace = 0.0, bottom = .11 )#.5
        plt.plot( x, data[:, 1:][ data[:, 0] == 2].T, '-x', color = 'gray' )
        plt.xlim( 0.5, 4.5 )
        plt.ylim( 0, 1.1 * max( data[:, 1:][ data[:, 0] == 2] ) )

        plt.xlabel( '$\mathrm{configuration}$', size = 20 )
        plt.xticks( ( 1, 2, 3, 4 ) ,
                     ( '$\mathrm{I}$', '$\mathrm{II}$', '$\mathrm{III}$', '$\mathrm{IV}$' ), size = 20, position = ( 0, -.01 ) )
        plt.ylabel( '$\mathrm{time\, [sec]}$', size = 20 )
        plt.yticks( ( 0, 50, 100, 150 ) ,
                     ( '$0$', '$50$', '$100$', '$150$' ) )
        plt.ylim( 0, 150 )
        plt.show()

    def plot_filament( self ):
        ''' Plot fillament_comb_int100_eps04_1.dat
        '''
        import matplotlib.pyplot as plt
        from matplotlib import rc
        rc( 'font', family = 'serif',
             style = 'normal', variant = 'normal', stretch = 'normal', size = 16 )

        data = self.data_class.data2.copy()
        data[data < 0] = self.data_class.data1[data < 0]
        x = [1, 2, 3, 4]
        plt.figure( 0 )
        plt.subplots_adjust( wspace = 0.0, hspace = 0.0, bottom = .11 )
        plt.plot( x, data[:, 1:][ data[:, 0] == 4].T, '-x', color = 'gray' )
        plt.xlim( 0.5, 4.5 )
        plt.ylim( 0, 1.1 * max( data[:, 1:][ data[:, 0] == 4] ) )

        plt.xlabel( '$\mathrm{configuration}$', size = 20 )
        plt.xticks( ( 1, 2, 3, 4 ) ,
                     ( '$\mathrm{I}$', '$\mathrm{II}$', '$\mathrm{III}$', '$\mathrm{IV}$' ), size = 20, position = ( 0, -.01 ) )
        plt.ylabel( '$\mathrm{time\, [sec]}$', size = 20 )
        plt.yticks( ( 0, 25, 50, 75, 100, 125 ) ,
                     ( '$0$', '$25$', '$50$', '$75$', '$100$', '$125$' ) )
        plt.ylim( 0, 125 )
        plt.annotate( '$E, A, \\theta, \lambda$',
                       xy = ( 3.4, 76 ), xycoords = 'data',
                       bbox = dict( boxstyle = "square,pad=.2", fc = "white" ),
                    xytext = ( 2.6, 85 ), #textcoords = 'offset points',
                    arrowprops = dict( arrowstyle = "->" ), size = 20,
                    )
        plt.annotate( '$E, A, \lambda, \\xi$\n$E, A, \\theta, \\xi$\n$ A, \\theta, \lambda, \\xi$\n$E, \\theta, \lambda, \\xi$',
                       xy = ( 3.45, 60 ), xycoords = 'data',
                       bbox = dict( boxstyle = "square,pad=.2", fc = "white" ),
                    xytext = ( 3.6, 15 ), #textcoords = 'offset points',
                    arrowprops = dict( arrowstyle = "->" ), size = 20,
                    )
        plt.show()

    def plot_filament_var2( self ):
        ''' Plot filament_comb_int5000_eps04_var2_1.dat
        '''
        import matplotlib.pyplot as plt
        from matplotlib import rc
        rc( 'font', family = 'serif',
             style = 'normal', variant = 'normal', stretch = 'normal', size = 16 )

        data = self.data_class.data2.copy()
        data[data < 0] = self.data_class.data1[data < 0]
        x = [1, 2, 3, 4]
        plt.figure( 0 )
        plt.subplots_adjust( wspace = 0.0, hspace = 0.0, bottom = .11 )
        plt.plot( x, data[:, 1:][ data[:, 0] == 2].T, '-x', color = 'gray' )
        plt.xlim( 0.5, 4.5 )
        plt.ylim( 0, 1.1 * max( data[:, 1:][ data[:, 0] == 2] ) )

        plt.xlabel( '$\mathrm{configuration}$', size = 20 )
        plt.xticks( ( 1, 2, 3, 4 ) ,
                     ( '$\mathrm{I}$', '$\mathrm{II}$', '$\mathrm{III}$', '$\mathrm{IV}$' ), size = 20, position = ( 0, -.01 ) )
        plt.ylabel( '$\mathrm{time\, [sec]}$', size = 20 )
        plt.yticks( ( 0, 50, 100, 150, 200 ) ,
                     ( '$0$', '$50$', '$100$', '$150$', '$200$' ) )
        plt.ylim( 0, 200 )
        plt.show()


    data_changed = Event( True )
    @on_trait_change( '+modified, data' )
    def _redraw( self ):
        figure = self.figure
        axes = figure.axes[0]
        axes.clear()
        x = [1, 2, 3, 4]

        if self.sub_slider_on == False:
            axes.plot( x, self.data[:, 1:][ self.data[:, 0] == self.comb_slider].T, '-x' )
        else:
            axes.plot( x, self.data[:, 1:][ self.data[:, 0] == self.comb_slider].T, '-x', color = 'gray' )
            axes.plot( x, self.data[:, 1:][ self.data[:, 0] == self.comb_slider][self.sub_slider, :].T, \
                       '-x', color = 'blue', linewidth = 2, label = self.data_text[ self.data[:, 0] == self.comb_slider][self.sub_slider][:-1] )
            axes.legend( loc = 0 )

        axes.set_xlim( 0.5, 4.5 )
        axes.set_ylim( 0, 1.1 * max( self.data[:, 1:][ self.data[:, 0] == self.comb_slider] ) )

        axes.set_xlabel( '\mathrm{configuration}' )
        axes.set_xticks( ( 1, 2, 3, 4 ) )
        axes.set_xticklabels( ( 'I', 'II', 'III', 'IV' ), position = ( 0, -.01 ) )
        self.data_changed = True

    view = View( Group( Item( 'data_class@', label = 'Files' ),
                        'data_sel',
                       Item( 'comb_slider' ),
                        HGroup( Item( 'sub_slider_on', label = 'sub_slider' ),
                                Item( 'sub_slider', show_label = False, springy = True, enabled_when = 'sub_slider_on == True' ) ),
                       Item( 'figure', style = 'custom',
                                  editor = MPLFigureEditor(),
                                  show_label = False ), ),
                resizable = True,
                scrollable = True,
                dock = 'tab',
                width = 0.8,
                height = 0.4 )





if __name__ == '__main__':
    data = Data()
    v = PView( data_class = data )

    v.configure_traits()

    #v.plot_pullout()
    #v.plot_pullout_var2()
    v.plot_filament()
    #v.plot_filament_var2()

