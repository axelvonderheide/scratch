
from enthought.traits.api import HasTraits, Directory, File
import os
from mfn_polar import MFnPolar

class LaTeXEngine( HasTraits ):
    '''Class managing output to a latex file.
    
    The interface of this class encapsulates the basic
    latex construct - title, author, section, table, table_row
    and image.
    '''
    
    #Initialization the variables
    image_plot_initial = 'polar_plot_initial.png'
    image_plot_fn1     = 'polar_fn1_plot.png'
    image_plot_fn2     = 'polar_fn2_plot.png'
    image_plot_uc      = 'polar_uc_plot.png'
    
    image_list =["image_plot_initial", "image_plot_fn1", "image_plot_fn2", "image_plot_uc"]
    
    title_page = 'Converting Polar Plot to PDF'
    author_name = 'Faeze'
    
    directory = Directory
    file = File

    def document_class(self, file):
        file.write('\documentclass[paper=a4, fleqn]{scrartcl}' + '\n')
    
    def use_package(self, file):
        file.write('\usepackage[latin1]{inputenc}' + '\n')
        file.write('\usepackage{dsfont}' + '\n')
        file.write('\usepackage{color}'+ '\n')
        file.write('\usepackage{amsmath}' + '\n')
        file.write('\usepackage{amssymb}' + '\n')
        file.write('\usepackage{comment}' + '\n')
        file.write('\usepackage{mathptmx}' + '\n')
        file.write('\usepackage[scaled=.90]{helvet}' + '\n')
        file.write('\usepackage{graphicx}' + '\n')
        file.write('\usepackage{eucal}' + '\n')
        file.write('\usepackage{palatino}' + '\n') 
        file.write('\usepackage{listings}' + '\n')   # Gives syntax highlighting for python code. 
        file.write('\usepackage{textcomp}' + '\n' + '\n' + '\n')

    def begin_document(self, file):
        file.write('\\begin{document}' + '\n' + '\n')    

    def end_document(self, file):
        file.write('\\end{document}' + '\n' )    
        
    def begin_title(self, file):
        file.write('\\begin{titlepage}' + '\n' + '\n')
        self.begin_center(file)
        file.write(self.title_page + '\n')
        self.end_center(file)
        
    def end_title(self, file):
        file.write('\\end{titlepage}' + '\n' + '\n')

    def author(self, file):
        file.write('\\author{' + self.author_name + '}' + '\n' + '\n')
        
    def section(self, file, name):
        file.write('\\section{' + name +'}' + '\n'  + '\n')
        
    def sub_section(self, file, name):
        file.write('\\subsection{' + name +'}' + '\n'  + '\n')
        
    def input_param(self, file, exported_obj ):
        print exported_obj.eparams
        eparams = exported_obj.eparams
        file.write('$')
        for e in eparams:
            file.write( e.replace("_","\_") + ' = ' + str( getattr( exported_obj, e ) )  + '\n' + '\\newline' + '\n')
        file.write('\\newline$' + '\n' + '\n') 
        
    def insert_image(self, file):
        image_list = self.image_list
        for l in image_list:
            name =  str( getattr( self, l ) )
            file.write('\\includegraphics[scale=0.6]{' + name + '}' + '\n' + '\n')
            
    def begin_center(self, file):
        file.write('\\begin{center}' + '\n')

    def end_center(self, file):
        file.write('\\end{center}' + '\n'  + '\n')
             
    def export_tex( self, exported_obj ):
        
        test_file = os.path.join('', 'Export_polar_latex')
        
        tex_file = test_file + '.tex'
        pdf_file = test_file + '.pdf'

        output_file = open( tex_file,'w')
        
        # get the parameters to be exported
        self.document_class( output_file )
        self.use_package(output_file)
        self.begin_document(output_file)
        self.author(output_file)
        self.begin_title(output_file)
        
        self.section(output_file, 'Parameters')
        self.sub_section(output_file, 'Initial Values')
        self.input_param(output_file, exported_obj )        
        
        self.section(output_file, 'show Plot')
        self.sub_section(output_file, 'Initial Plot, Function Polar plots and Unit Circle Plot')
        self.insert_image(output_file)
        self.end_title(output_file)
        self.end_document(output_file)
        
        output_file.close()
        
        os.system( 'pdflatex ' + tex_file )
        os.system( 'okular ' + pdf_file )

lengine = LaTeXEngine()
# setup of the output

pp = MFnPolar()
#pp.export_latex( "file" )

lengine.export_tex( pp )

#LaTeXEngine().export_tex( MFnPolar() )

