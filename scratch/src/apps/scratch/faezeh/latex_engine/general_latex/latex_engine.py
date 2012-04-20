
from enthought.traits.api import HasTraits, Directory, File
import os

class LaTeXEngine( HasTraits ):
    '''Class managing output to a latex file.
    
    The interface of this class encapsulates the basic
    latex construct - title, author, section, table, table_row
    and image.
    '''
    
    #Initialization the variables
    height = 10.0
    width  = 15.0

    eparams = ["height", "width" ]

    image1 = 'noninteractive.png'
    
    image_list =["image1"]
    
    title_page = 'Converting to PDF'
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
        
    def input_param(self, file):
        eparams = self.eparams
        for e in eparams:
            file.write(e + ' = $' + str( getattr( self, e ) ) +'$' + '\n' + '\n')

    def insert_image(self, file):
        image_list = self.image_list
        for l in image_list:
            name =  str( getattr( self, l ) )
            file.write('\\includegraphics[scale=0.6]{' + name + '}' + '\n' + '\n')
            
    def begin_center(self, file):
        file.write('\\begin{center}' + '\n')

    def end_center(self, file):
        file.write('\\end{center}' + '\n'  + '\n')
             
def export_tex( exported_obj ):
    # get the parameters to be exported
    test_file = os.path.join('', 'Export_file_latex')
    
    output_file = open(test_file + '.tex','w')
    exported_obj.document_class( output_file )
    exported_obj.use_package(output_file)
    exported_obj.begin_document(output_file)
    exported_obj.author(output_file)
    exported_obj.begin_title(output_file)
    
    exported_obj.section(output_file, 'Parameters')
    exported_obj.input_param(output_file)
    exported_obj.section(output_file, 'show Plot')
    exported_obj.insert_image(output_file)
    exported_obj.end_title(output_file)
    exported_obj.end_document(output_file)
    
    output_file.close()

ip = LaTeXEngine()

export_tex( ip )
   
