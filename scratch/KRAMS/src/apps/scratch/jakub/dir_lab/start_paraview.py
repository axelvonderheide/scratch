'''
Created on May 11, 2009

@author: jakub
'''
    start_paraview = Button('Paraview')
    
    @on_trait_change('start_paraview')
    def start_pw(self,event = None):
        os.system("paraview --data="+self.var+"%(direction)d%(pos)s_..vtk" \
                  %{'direction': (self.idx+1), "pos": self.position})