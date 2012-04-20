'''
Created on Apr 9, 2010

@author: alexander
'''

from os.path import \
    join
    
from string import strip
import os

import csv

from numpy import array

from scipy.io import read_array
from numpy import loadtxt, argmax, polyfit, poly1d, frompyfunc


from enthought.traits.ui.file_dialog  \
    import open_file, FileInfo, TextInfo, ImageInfo

import sys

from matplotlib.mlab import csv2rec, csvformat_factory


from numpy import *

from scipy import loadtxt

from time import *

from os.path import join

from enthought.traits.ui.api import \
    View, Item, FileEditor, HSplit, Group, VSplit, CancelButton, OKButton, \
    Handler

from promod.exdb.ex_run import ExRun
from promod.exdb.ex_run_view import ExRunView
data_file_editor = FileEditor( filter = ['*.DAT'] )

from promod.simdb import SimDB
simdb = SimDB()

#test_file = join( simdb.exdata_dir, 
#                  'tensile_tests', 
#                  'TT-10a',
#                  'TT11-10a-V1.DAT' )
#
#ex_run = ExRun( test_file )
#
#print 'ex_run', ex_run.ex_type.E_c
#print 'ex_run', ex_run.


file_name = '/home/alexander/simdb/exdata/plate_tests/PT-7u/PT07-7u.ASC'
#file_name = '/home/alexander/simdb/exdata/plate_tests/PT-10a/PT11-10a_original.ASC'




data_arr = loadtxt( file_name, 
                    delimiter = ';')

print 'data_arr', data_arr

def loadtxt_novalue( file_name ):
    '''Return an data array similar to loadtxt. 
    "NOVALUE" entries are replaced by the value of the previous line.
    '''
    file = open( file_name,'r')
    lines = file.readlines()
    n_columns = len( lines[0].split(';') )
    data_array = zeros( n_columns )
    n = 0
    for line in lines:
        line_split =  lines[n].split(';') 
        m = 0 
        for value in line_split:
            if value == 'NOVALUE':
                print 'NOVALUE entry in line', n, 'position', m, 'found and replaced'
                line_split[m] = lines[n-1].split(';')[m] 
            m += 1
        line_array = array(line_split, dtype = float)
        data_array = vstack([ data_array, line_array ])
        n += 1
    return data_array




def check_for_novalue_entries( file_name ):
    '''check if the data file contains "NOVALUE" entries.
    Returns the total number of NOVALUE entries and print their position.
    '''
    n_novalue_entries = 0
    file = open( file_name,'r')
    lines = file.readlines()
    n_columns = len( lines[0].split(';') )
    n = 0
    for line in lines:
        line_split = lines[n].split(';') 
        m = 0 
        for value in line_split:
            if value == 'NOVALUE':
                print 'NOVALUE entry in line', n, 'position', m, 'found'
                n_novalue_entries += 1
            m += 1
        n += 1
    print 'total number of NOVALUE entries =', n_novalue_entries
    return n_novalue_entries

t1 = time()
n_novalue_entries = check_for_novalue_entries( file_name )
t2 = time()
Dt = t2-t1
print 'Dt', Dt, '\n'


t1 = time()
cc = loadtxt_novalue( file_name )
t2 = time()
Dt = t2-t1
print 'Dt',Dt, '\n'

print 'cc' , cc








#"".join(['a','b'])


##aaa = csv2rec( file_name )
##print 'aaa', aaa.shape
##print 'aaa', array( aaa[0][0].split(';') )
#
#file = open( file_name,'r')
#lines = file.readlines()
#print 'aaa', lines#.split(';')



#.sub


#
#def Xloadtxt_novalue( file_name ):
#    '''Return an data array similar to loadtxt. 
#    "NOVALUE" entries are replaced by the value of the previous line.
#    '''
#    file = open( file_name,'r')
#    lines = file.readlines()
#    n_columns = len( lines[0].split(';') )
#    data_array = zeros( n_columns )
#    n = 0
#    for line in lines:
#        line_split =  lines[n].split(';') 
#        try:
#            line_array = array(line_split, dtype = float)
#        except ValueError:
#            m = 0 
#            for value in line_split:
#                if value == 'NOVALUE':
#                    print 'NOVALUE entry in line', n, 'position', m, 'found and replaced'
#                    line_split[m] = lines[n-1].split(';')[m] 
#                m += 1
#            line_array = array(line_split, dtype = float)
#        data_array = vstack([ data_array, line_array ])
#        n += 1
#    return data_array
#
#t1 = time()
#cc = Xloadtxt_novalue( file_name )
#t2 = time()
#Dt = t2-t1
#print 'Dt',Dt

#lines = file.read(14).split()
#print 'lines', len(lines[0].split(';'))
#for n in [1]:
#print 'll',ll
#
#
#
#
#data_array = loadtxt( file_name , 
#                      delimiter = ';' )
