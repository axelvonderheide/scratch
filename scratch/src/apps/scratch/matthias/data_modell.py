'''
Created on June 09, 2011

@author: schmerl
'''
from enthought.traits.api import HasTraits, Property, cached_property, Event, \
    Array, Instance, Int, Directory, Range, on_trait_change, Bool, Trait, Constant, \
    Tuple, Interface, implements, Enum
from enthought.traits.trait_types import DelegatesTo
from enthought.traits.ui.api import Item, View, HGroup, RangeEditor
from numpy import loadtxt, min, array, arange, ones_like, cumsum, vstack, \
    hstack, sum, zeros_like, zeros, ones, where, unique, pi, invert, \
    prod
    
import math




class DataModell( HasTraits ):
    '''
        Data_Modell of the creasepattern
    '''
    
    # just 3 Arrays representing the Vertice
    # TODO BETTER!!
    # take one Array of 3DPoints 
    # Crease lines has to be defined with mountain/valley information
    # Maybe creating the mesh already in here?
    
    
    x = array([0,-1,2,0,1,-1,1,2])
    y = array([0,2,2,1,1,-1,0,-1])
    z = array([0,1,1,0,0,1,0,1])

    

class Crease_line( HasTraits ):
    '''
        Datas of a crease-line
        described by 2 vertices
            assuming the vertices describing an unfolded creasepattern they 
            are in the x-y plane, so z will be in all vertives the same
        and the crease angel roh
            roh is an angle between -pi<roh<pi
            roh = zero for non folded paper
        alpha is the rotation angel in the x-y Plane, and can be calculated 
        by the vertices (or otherway round)
            0<alpha<2PI
            alpha is always messured from the x-achses
    '''
    #Initialise
    vertix1 =array([0,0,0])
    vertix2 =array([0,0,0])
    alpha = 0.0
    roh = 0.0
    
    def calc_alpha(self):
        deltaY = self.vertix2[1]-self.vertix1[1]   
        deltaX = self.vertix2[0]-self.vertix1[0]
        self.alpha = math.atan(deltaY/deltaX)
        while(self.alpha > 2*math.pi):
            self.alpha -= 2*math.pi
            
        while(self.alpha < 0):
            self.alpha += 2*math.pi
    
    def calc_vertix2(self):
        y = math.tan(self.alpha)
        self.vertix2 = array([self.vertix1[0]+1,self.vertix1[1]+y,self.vertix1[2]])
            
    def get_vertix(self,x,y,z):
        self.vertix1= array([x,y,z])
        self.calc_alpha()
        
    def get_neighbor(self,x,y,z):
        self.vertix2 = array([x,y,z])
        self.calc_alpha()
        
    def get_alpha(self,a):
        while(a >2*math.pi):
            a -= 2*math.pi
        while(a<0):
            a += 2*math.pi
            
        self.alpha = a
        self.calc_vertix2()
        
