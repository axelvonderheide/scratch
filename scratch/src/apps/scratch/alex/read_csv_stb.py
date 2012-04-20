'''
Created on Jun 28, 2010

@author: alexander
'''

from numpy import array, loadtxt, arange, sqrt, zeros, arctan, sin, cos, ones_like, \
                vstack, savetxt, hstack, argsort, fromstring, zeros_like, \
                copy, c_, newaxis, argmax, where, argsort, sqrt, frompyfunc

file_name = 'Schale_Knotenkoordinaten_orig.csv'

def read_node_coords_stb( file_name ):
    '''read the nodal coordinates of the mid-surface
    defined in a csv-file. To export the excel sheet 
    to csv use ";" as a field delimiter and "" (none) 
    as a text delimiter. 
    Note that some lines do not contain values!
    '''

    file = open( file_name, 'r' )

    # read the column headings (first two lines)
    #
    first_line = file.readline()
    second_line = file.readline()
    column_headings = second_line.split( ';' )
    # remove '\n' from last string element in list
    column_headings[-1] = column_headings[-1][:-1]
    column_headings_arr = array( column_headings )

    # check in which column the node number and the 
    # carthesian coordinates can be found 
    #
    elem_no_idx = where( 'Nr.' == column_headings_arr )[0]
    X_idx = where( 'X [m]' == column_headings_arr )[0]
    Y_idx = where( 'Y [m]' == column_headings_arr )[0]
    Z_idx = where( 'Z [m]' == column_headings_arr )[0]

    lines = file.readlines()

    lines_list = [ line.split( ';' ) for line in lines ]

    ll = []
    for line in lines_list:
        # check if lien contains values or only a node number!
        #
        if line[1] == 'Standard':
            ll.append( [line[elem_no_idx], line[X_idx], line[Y_idx], line[Z_idx]] )

    input_arr = array( ll, dtype = 'float_' )

    node_no = input_arr[:, 0]
    X = input_arr[:, 1]
    Y = input_arr[:, 2]
    Z = input_arr[:, 2]

    return node_no, X, Y, Z

node_no, X, Y, Z = read_node_coords_stb( file_name )

print c_[ node_no, X, Y, Z ].shape



