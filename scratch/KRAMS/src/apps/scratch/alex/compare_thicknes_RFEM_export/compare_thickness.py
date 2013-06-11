'''
Created on Jun 16, 2010

@author: andreas
'''

from enthought.traits.api import \
    HasTraits, Float, Array, implements, Property, cached_property, Instance, Enum, \
    Dict, Bool, Int

from numpy import \
    array, tensordot, dot, zeros, c_, ix_, mgrid, arange, \
    where, sum, sin, cos, vstack, hstack, argmax, newaxis, size, \
    shape, sqrt, frompyfunc, ones_like, loadtxt, arange, sqrt, \
    zeros, arctan, sin, cos, ones_like, \
    vstack, savetxt, hstack, argsort, fromstring, zeros_like, \
    copy, c_, where, mean, arctan, cos, min, argmin

from os.path import join

from math import pi

from enthought.traits.api import \
    HasTraits, Float, Array, Bool, Enum, Dict

from numpy import \
    c_, ix_, mgrid, transpose, shape

from rsurface_reader import \
    read_rsurface, normalize_rsurfaces

# Interpolation
from scipy.interpolate import Rbf

from math import pi

import csv

from enthought.mayavi.mlab import colorbar, show, points3d
from enthought.mayavi.api import Engine

from enthought.mayavi import \
    mlab

compare = True

file_name = 'RFEM_file_export/output_data_thickness.csv'
file_name_wf = 'RFEM_file_export/output_data_thickness_without_filter.csv'

input_arr = loadtxt( file_name, delimiter = ';', skiprows = 1 )
input_arr_wf = loadtxt( file_name_wf, delimiter = ';', skiprows = 1 )

#node_number;x[m];y[m];z[m]

n_elem = input_arr[:, 0]
X = input_arr[:, 1]
Y = input_arr[:, 2]
t = input_arr[:, 3]

n_elem_wf = input_arr_wf[:, 0]
X_wf = input_arr_wf[:, 1]
Y_wf = input_arr_wf[:, 2]
t_wf = input_arr_wf[:, 3]

error_X_abs = X_wf - X
error_Y_abs = Y_wf - Y
error_t_abs = abs( t_wf - t )

print 'max(error_X_abs)', max( error_X_abs )
print 'error_abs', error_X_abs

print 'max(error_Y_abs)', max( error_Y_abs )
print 'error_abs', error_Y_abs

print 'max(error_t_abs)', max( error_t_abs )
print 'error_abs', error_t_abs
print 'X,Y pos max(error_t_abs)', X[ argmax( error_t_abs ) ], Y[ argmax( error_t_abs ) ]

error_t_rel = abs( 1 - t / t_wf ) * 100

print 'max(error_t_rel)', max( error_t_rel )
print 'X,Y pos max(error_t_rel)', X[ argmax( error_t_rel ) ], Y[ argmax( error_t_rel ) ]
print 'error_t_rel', error_t_rel

Z = zeros_like( X )

@show
def plot_XYZ():
    mlab.points3d( X, Y, Z,
                   error_t_rel,
#                   error_t_abs,
                   colormap = "YlOrBr",
                   mode = "cube",
                   scale_factor = 0.015 )
    mlab.scalarbar( orientation = 'vertical' )
#
#plot_XYZ()

#--------------------------
# mid surface:
#--------------------------

file_name = 'RFEM_file_export/output_data_midsurface.csv'
file_name_wf = 'RFEM_file_export/output_data_midsurface_without_filter.csv'

input_arr = loadtxt( file_name, delimiter = ';', skiprows = 1 )
input_arr_wf = loadtxt( file_name_wf, delimiter = ';', skiprows = 1 )

#node_number;x[m];y[m];z[m]

n_node = input_arr[:, 0]
X = input_arr[:, 1]
Y = input_arr[:, 2]
Z = input_arr[:, 3]

n_node_wf = input_arr_wf[:, 0]
X_wf = input_arr_wf[:, 1]
Y_wf = input_arr_wf[:, 2]
Z_wf = input_arr_wf[:, 3]

error_X_abs = X_wf - X
error_Y_abs = Y_wf - Y
error_Z_abs = abs( Z_wf - Z )

print 'max(error_X_abs)', max( error_X_abs )
print 'error_abs', error_X_abs

print 'max(error_Y_abs)', max( error_Y_abs )
print 'error_abs', error_Y_abs

print 'max(error_Z_abs)', max( error_Z_abs )
print 'error_abs', error_Z_abs
print 'X,Y pos max(error_Z_abs)', X[ argmax( error_Z_abs ) ], Y[ argmax( error_Z_abs ) ]

error_Z_rel = where( Z_wf != 0, abs( 1 - Z / Z_wf ) * 100, 0 )

print 'max(error_Z_rel)', max( error_Z_rel )
print 'X,Y pos max(error_Z_rel)', X[ argmax( error_Z_rel ) ], Y[ argmax( error_Z_rel ) ]
print 'error_Z_rel', error_Z_rel

@show
def plot_XYZ():
    mlab.points3d( X, Y, Z_wf,
                   error_Z_abs,
                   colormap = "YlOrBr",
                   mode = "cube",
                   scale_factor = 100. )
    mlab.scalarbar( orientation = 'vertical' )

plot_XYZ()















# find the lines in the midsurface file, that do not contain values
#
#file_name = 'input_data_midsurface.csv'
#file = open( file_name, 'r' )
#
## read first two lines
##
#first_line = file.readline()
#second_line = file.readline()
#
## read remaining lines
##
#lines = file.readlines()
#
#lines_list = [ line.split( ';' ) for line in lines ]
#
#empty_idx_list = []
#for idx, line in enumerate( lines_list ):
#    # check if lien contains values or only a node number!
#    #
#    print 'line[1]', line[1]
#    if line[1] != 'Standard':
#        empty_idx_list.append( idx )
#
#print 'empty_idx_list', empty_idx_list
















