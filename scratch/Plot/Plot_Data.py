from etsproxy.traits.api import \
    HasTraits, Instance, Int, Array, List, \
    cached_property, Property, Float, Bool, String
from operator import attrgetter
import numpy as np

import copy
from spirrid.rv import RV
from matplotlib import pyplot as plt
import time as t
import mayavi
import numpy
import pickle
import os
class Plot_view( HasTraits ):
        
        def plot_control( self, pc, ps, plot_es, plot_sw, plot_hist, plot_hist_sigmas, plot_folder, form_es, form_sw ):
            if plot_es:
                plt.figure()
                lw1, ls1, c1, lw2, ls2, c2 = form_es
                for foldername in plot_folder:
                    self.open_folder( foldername )
                    self.plot_e_s_f( pc, ps, lw1, ls1, c1, lw2, ls2, c2 )
                    self.close_folder()
                if pc:plt.plot( 0, 0, linewidth = lw1, linestyle = ls1, color = c1, label = 'combined' )
                if ps:plt.plot( 0, 0, linewidth = lw2, linestyle = ls2, color = c2, label = 'solo' )
                plt.xlabel( 'composite strain [-]' )
                plt.ylabel( 'composite stress [MPa]' )
                plt.legend( loc = 'best' )
                    
            if plot_sw:
                plt.figure()
                lw1, ls1, c1, lw2, ls2, c2, meanplot, stdevplot, maxplot = form_sw
                for foldername in plot_folder:
                    self.open_folder( foldername )
                    self.plot_s_w_f( pc, ps, lw1, ls1, c1, lw2, ls2, c2, meanplot, stdevplot, maxplot, )
                    self.close_folder()   
                if pc:plt.plot( 0, 0, linewidth = lw1, linestyle = ls1, color = c1, label = 'combined' )
                if ps:plt.plot( 0, 0, linewidth = lw2, linestyle = ls2, color = c2, label = 'solo' )
                plt.xlabel( 'composite stress [MPa]' )
                plt.ylabel( 'mean crack opening [mm]' )
                    
            if plot_hist:
                plt.figure()
                for foldername in plot_folder:
                    self.open_folder( foldername )
                    self.plot_hist( pc, ps, plot_hist_sigmas )
                    self.close_folder()
                plt.xlabel( 'crack opening [mm]' )
                plt.ylabel( 'number of cracks' )
                plt.legend( loc = 'best' )
                
            plt.show()
                    
        def open_folder( self, folder ):
                name = 'Data/{}'.format( folder )
                os.chdir( name )
        def close_folder( self ):
                os.chdir( os.pardir )
                os.chdir( os.pardir )
            
        def plot_e_s_f( self, pc, ps, lw1, ls1, c1, lw2, ls2, c2 ):
                if pc:
                    combined_es = open( 'combined_es.pkl', 'rb' )
                    combined_e_s_arr = pickle.load( combined_es )
                    combined_es.close()
                    plt.plot( combined_e_s_arr[0], combined_e_s_arr[1], linewidth = lw1, linestyle = ls1, \
                               color = c1 )
                if ps:
                    solo_es = open( 'solo_es.pkl', 'rb' )
                    solo_e_s_arr = pickle.load( solo_es )
                    solo_es.close()
                    plt.plot( solo_e_s_arr[0], solo_e_s_arr[1], linewidth = lw2, \
                              linestyle = ls2, color = c2 ) 

        def plot_s_w_f( self, pc, ps, lw1, ls1, c1, lw2, ls2, c2, meanplot, stdevplot, maxplot ):
                if pc:
                    combined_sw = open( 'combined_sw.pkl', 'rb' )
                    combined_sw_arr = pickle.load( combined_sw )
                    combined_sw.close()
                    if meanplot:plt.plot( combined_sw_arr[0], combined_sw_arr[1] , linewidth = lw1, linestyle = ls1, \
                               color = c1 )
                    if stdevplot:
                        plt.plot( combined_sw_arr[0], combined_sw_arr[1] - combined_sw_arr[3] , linewidth = lw1 / 2., linestyle = '--', \
                               color = c1 )
                        plt.plot( combined_sw_arr[0], combined_sw_arr[1] + combined_sw_arr[3] , linewidth = lw1 / 2., linestyle = '--', \
                               color = c1 )
                    if maxplot:plt.plot( combined_sw_arr[0], combined_sw_arr[2] , linewidth = lw1 / 2., linestyle = '--', \
                               color = 'r' )
                if ps:
                    solo_sw = open( 'solo_sw.pkl', 'rb' )
                    solo_sw_arr = pickle.load( solo_sw )
                    solo_sw.close()
                    if meanplot:plt.plot( solo_sw_arr[0], solo_sw_arr[1], linewidth = lw2, \
                              linestyle = ls2, color = c2 ) 
                    if stdevplot:
                        plt.plot( solo_sw_arr [0], solo_sw_arr [1] - solo_sw_arr[3] , linewidth = lw1 / 2., linestyle = '--', \
                               color = c1 )
                        plt.plot( solo_sw_arr[0], solo_sw_arr[1] + solo_sw_arr[3] , linewidth = lw1 / 2., linestyle = '--', \
                               color = c1 )
                    if maxplot:plt.plot( solo_sw_arr[0], solo_sw_arr[2] , linewidth = lw1 / 2., linestyle = '--', \
                               color = 'r' )
                
        def plot_hist( self, pc, ps, plot_hist_sigmas ):
            if pc:
                combined_hist = open( 'combined_hist.pkl', 'rb' )
                combined_hist_list = pickle.load( combined_hist )
                combined_hist.close()
                combined_es = open( 'combined_es.pkl', 'rb' )
                combined_e_s_arr = pickle.load( combined_es )
                combined_es.close()
                combined_e_s_arr
                for sigma in plot_hist_sigmas:
                    idx = np.argmin( np.abs( combined_e_s_arr[1] - sigma ) )
                    plt.hist( combined_hist_list[idx], bins = 20, label = 'load = {} MPa'.format( sigma ) )
            
if __name__ == '__main__':
    # Control
    combined = True
    solo = False
    plot_es = True
    
    plot_sw = False
    meanplot = True
    stdevplot = True
    maxplot = True
    plot_hist = False
    plot_hist_sigmas = [ 5., 7., 8.]
    plot_folder = [ '0.005', '0.015']
    ######
    form_es = [2., '-', 'k', 2., '-', 'b']
    form_sw = [2., '-', 'k', 2., '-', 'b', meanplot, stdevplot, maxplot]
    ini = Plot_view()
    ini.plot_control( combined, solo, plot_es, plot_sw, plot_hist, plot_hist_sigmas, plot_folder, form_es, form_sw )

                
