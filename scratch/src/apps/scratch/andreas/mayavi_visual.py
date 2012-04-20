'''
Created on 07.07.2010

@author: Andy
'''
    plot_column = Enum( values = 'columns' )
    plot = Button
    def _plot_fired( self ):
        X = self.info_shell.X[:, 0]
        Y = self.info_shell.Y[:, 0]
        Z = self.info_shell.Z[:, 0]
        plot_col = getattr( self, self.plot_column )[:, 0]
        if self.plot_column == 'n_tex':
            plot_col = where(plot_col < 0, 0, plot_col)
            
        print '*** plotting data***'
        mlab.figure( figure = "Demonstrator", bgcolor=(1.0,1.0,1.0),fgcolor=(0.0,0.0,0.0) )
        mlab.points3d( X, Y, (-1.0)*Z, plot_col, colormap = "gist_rainbow", mode = "cube", scale_factor = 0.1 )
        mlab.scalarbar( title= self.plot_column, orientation = 'vertical')