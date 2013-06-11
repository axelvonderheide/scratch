from numpy import array, loadtxt, arange, sqrt, zeros, arctan, sin, cos, ones_like,\
                vstack, savetxt, hstack, argsort, fromstring, transpose
from math import pi
from string import split
from input_params import D_b, D_t, f_ctk, f_m


# ------------------------------------------------------------
# Read input txt-file and get data into numpy format
# ------------------------------------------------------------
# open input file containing raw data:
file = open('Rohdaten_GdG','r')

# open subsidary file for the conversion between "." and "," defined floats:
a_data = open('adapted_data','w')

# search and replace "," with ".":
adapted = file.read().replace(",",".")

# print 'adapted':
a_data.write(adapted)
a_data.close()

# get subsidary file to retriev input data as array of floats:
input_arr = loadtxt( 'adapted_data' , usecols=(-8,-7,-6,-5,-4,-3), delimiter='\t' )

file = open( 'adapted_data','r')

# moments
mx  = input_arr[:,0]
my  = input_arr[:,1]
mxy = input_arr[:,2]


# normalforces
nx  = input_arr[:,3]
ny  = input_arr[:,4]
nxy = input_arr[:,5]


# ------------------------------------------------------------
# Get element number and node number
# ------------------------------------------------------------
elem_no = ones_like(mx)
node_no  = ones_like(mx)

cr_elem = 0
for i in range( elem_no.shape[0] ):
    line = file.readline()
    if line.isspace():
        line = file.readline()
    value = split( line, '\t' )
    if value[0]:
        cr_elem = int(value[0])
        cr_node = int(value[1] )
    elem_no[i] = cr_elem
    node_no[i] = cr_node
    #if not line: break
print "elem_no ",elem_no
    
row_no = arange(elem_no.shape[0])

# ------------------------------------------------------------
# Get area A, moment of inertia W
# ------------------------------------------------------------
# derive the element thickness based on the element number
# (the division by 100 (=int) yields an integer)
D_elem = D_b - (elem_no/100-1)/11 * (D_b - D_t)

A = D_elem * 1.
W = 1. * D_elem**2/6.


# ------------------------------------------------------------
# Index M: principle moments with corresponding normalforces
# ------------------------------------------------------------

# principle moments
m1 = 0.5*(mx+my) + 0.5*sqrt((mx-my)**2 + 4*mxy**2)
m2 = 0.5*(mx+my) - 0.5*sqrt((mx-my)**2 + 4*mxy**2)

# principle angle of the moments:
alpha_M = pi/2*ones_like(m1)
alpha_M[m1 != m2] = arctan(mxy/(m2-mx))

# transform the values in the principle direction:
mx_M = 0.5*(mx+my) - 0.5*(mx-my)*cos(2*alpha_M) - mxy*sin(2*alpha_M)
my_M = 0.5*(mx+my) + 0.5*(mx-my)*cos(2*alpha_M) + mxy*sin(2*alpha_M)

nx_M = 0.5*(nx+ny) - 0.5*(nx-ny)*cos(2*alpha_M) - nxy*sin(2*alpha_M)
ny_M = 0.5*(nx+ny) + 0.5*(nx-ny)*cos(2*alpha_M) + nxy*sin(2*alpha_M)


# ------------------------------------------------------------
# Index N: principle normal forces with corresponding moments
# ------------------------------------------------------------

# principle normal forces
n1 = 0.5*(nx+ny) + 0.5*sqrt((nx-ny)**2 + 4*nxy**2)
n2 = 0.5*(nx+ny) - 0.5*sqrt((nx-ny)**2 + 4*nxy**2)

# principle angle of the moments:
alpha_N = pi/2*ones_like(n1)
alpha_N[n1 != n2] = arctan(nxy/(n2-nx))

nx_N = 0.5*(nx+ny) - 0.5*(nx-ny)*cos(2*alpha_N) - nxy*sin(2*alpha_N)
ny_N = 0.5*(nx+ny) + 0.5*(nx-ny)*cos(2*alpha_N) + nxy*sin(2*alpha_N)

mx_N = 0.5*(mx+my) - 0.5*(mx-my)*cos(2*alpha_N) - mxy*sin(2*alpha_N)
my_N = 0.5*(mx+my) + 0.5*(mx-my)*cos(2*alpha_N) + mxy*sin(2*alpha_N)


# ------------------------------------------------------------
# Evaluation of stresses
# ------------------------------------------------------------ 
# Index X: evaluation of stresses in transformed X-direction
# Index Y: evaluation of stresses in transformed Y-direction

# case: MX
sig_n_MX = nx_M / A
sig_m_MX = mx_M / W
eta_n_MX = sig_n_MX / f_ctk
eta_m_MX = sig_m_MX / f_m
eta_tot_MX = eta_n_MX + eta_m_MX
MX_ix = argsort( eta_tot_MX )

# case: MY
sig_n_MY = nx_M / A
sig_m_MY = mx_M / W
eta_n_MY = sig_n_MY / f_ctk
eta_m_MY = sig_m_MY / f_m
eta_tot_MY = eta_n_MY + eta_m_MY
MY_ix = argsort( eta_tot_MY )

# case: NX
sig_n_NX = nx_N / A
sig_m_NX = mx_N / W
eta_n_NX = sig_n_NX / f_ctk
eta_m_NX = sig_m_NX / f_m
eta_tot_NX = eta_n_NX + eta_m_NX
NX_ix = argsort( eta_tot_NX )

# case: NY
sig_n_NY = nx_N / A
sig_m_NY = mx_N / W
eta_n_NY = sig_n_NY / f_ctk
eta_m_NY = sig_m_NY / f_m
eta_tot_NY = eta_n_NY + eta_m_NY
NY_ix = argsort( eta_tot_NY )

max_ix = 4


# ------------------------------------------------------------
# MX - Save the output-file:
# ------------------------------------------------------------ 

max_ix = 4

save_arr_MX = transpose(vstack(( \
                               row_no[ MX_ix[:max_ix]],  \
                               elem_no[ MX_ix[:max_ix]], \
                               node_no[ MX_ix[:max_ix]], \
                               mx[ MY_ix[:max_ix]], \
                               my[ MY_ix[:max_ix]], \
                               mxy[ MY_ix[:max_ix]], \
                               nx[ MY_ix[:max_ix]], \
                               ny[ MY_ix[:max_ix]], \
                               nxy[ MY_ix[:max_ix]], \
                               m1[ MX_ix[:max_ix]], \
                               m2[ MX_ix[:max_ix]], \
                               alpha_M[ MX_ix[:max_ix]], \
                               mx_M[ MX_ix[:max_ix]],\
                               my_M[ MX_ix[:max_ix]], \
                               nx_M[ MX_ix[:max_ix]], \
                               ny_M[ MX_ix[:max_ix]], \
                               sig_n_MX[ MX_ix[:max_ix]], \
                               sig_m_MX[ MX_ix[:max_ix]], \
                               eta_n_MX[ MX_ix[:max_ix]], \
                               eta_m_MX[ MX_ix[:max_ix]], \
                               eta_tot_MX[ MX_ix[:max_ix]]
                               )))
                                        
save_list_MX = save_arr_MX.tolist()

# specify the labels of the outputfile:
label_list_MX = ['row_no', \
              'elem_no', \
              'node_no', \
              'mx', \
              'my', \
              'mxy', \
              'nx' ,\
              'ny', \
              'nxy', \
              'm1', \
              'm2', \
              'alpha_M', \
              'mx_M', \
              'my_M', \
              'nx_M', \
              'ny_M', \
              'sig_n_MX', \
              'sig_m_MX', \
              'eta_n_MX', \
              'eta_m_MX', \
              'eta_tot_MX'
              ]

# save the array input to file 'report' applying formated strings
output_GdG_file = open('output_GdG','w')

caption_str_MX = str('GdG : Evaluation for case MX:')
label_str_MX = str()
dash_line_str = '\n'
for i in range(len(label_list_MX)):
    label_str_MX = label_str_MX + '%010s '
    dash_line_str = dash_line_str + '---------- '
    
header_str_MX = caption_str_MX + dash_line_str + '\n' + label_str_MX +  dash_line_str + '\n'

output_GdG_file.write(header_str_MX %tuple(label_list_MX))

for line in save_list_MX:
    output_GdG_file.write('%10.0f %10.0f %10.0f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f \n' \
                     %(line[0], line[1], line[2], line[3], line[4], line[5], line[6], line[7], line[8], line[9], line[10], line[11], line[12], line[13], line[14], line[15], line[16], line[17], line[18], line[19], line[20] ) )

output_GdG_file.close()

# ------------------------------------------------------------
# MY - Save the output-file:
# ------------------------------------------------------------ 

max_ix = 4

save_arr_MY = transpose(vstack(( row_no[ MY_ix[:max_ix]], \
                       elem_no[ MY_ix[:max_ix]], \
                       node_no[ MY_ix[:max_ix]],\
                       mx[ MY_ix[:max_ix]], \
                       my[ MY_ix[:max_ix]], \
                       mxy[ MY_ix[:max_ix]], \
                       nx[ MY_ix[:max_ix]], \
                       ny[ MY_ix[:max_ix]], \
                       nxy[ MY_ix[:max_ix]], \
                       ###
                       m1[ MY_ix[:max_ix]], \
                       m2[ MY_ix[:max_ix]], \
                       alpha_M[ MY_ix[:max_ix]], \
                       mx_M[ MY_ix[:max_ix]],\
                       my_M[ MY_ix[:max_ix]], \
                       nx_M[ MY_ix[:max_ix]], \
                       ny_M[ MY_ix[:max_ix]], \
                       # MY
                       sig_n_MY[ MY_ix[:max_ix]], \
                       sig_m_MY[ MY_ix[:max_ix]], \
                       eta_n_MY[ MY_ix[:max_ix]], \
                       eta_m_MY[ MY_ix[:max_ix]], \
                       eta_tot_MY[ MY_ix[:max_ix]] )))
                                        
save_list_MY = save_arr_MY.tolist()

# specify the labels of the outputfile:
label_list_MY = ['row_no', \
              'elem_no', \
              'node_no', \
              'mx', \
              'my', \
              'mxy', \
              'nx' ,\
              'ny', \
              'nxy', \
              'm1', \
              'm2', \
              'alpha_M', \
              'mx_M', \
              'my_M', \
              'nx_M', \
              'ny_M', \
              'sig_n_MY', \
              'sig_m_MY', \
              'eta_n_MY', \
              'eta_m_MY', \
              'eta_tot_MY'
              ]

# save the array input to file 'report' applying formated strings
output_GdG_file = open('output_GdG','a')

caption_str_MY = str('GdG : Evaluation for case MY:')
label_str_MY = str()
dash_line_str = '\n'
for i in range(len(label_list_MY)):
    label_str_MY = label_str_MY + '%010s '
    dash_line_str = dash_line_str + '---------- '
    
header_str_MY = '\n \n \n' + caption_str_MY + dash_line_str + '\n' + label_str_MY +  dash_line_str + '\n'

output_GdG_file.write(header_str_MY %tuple(label_list_MY))

for line in save_list_MY:
    output_GdG_file.write('%10.0f %10.0f %10.0f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f \n' \
                     %(line[0], line[1], line[2], line[3], line[4], line[5], line[6], line[7], line[8], line[9], line[10], line[11], line[12], line[13], line[14], line[15], line[16], line[17], line[18], line[19], line[20] ) )

output_GdG_file.close()



# ------------------------------------------------------------
# NX - Save the output-file:
# ------------------------------------------------------------ 

max_ix = 4

save_arr_NX = transpose(vstack(( row_no[ NX_ix[:max_ix]], \
                       elem_no[ NX_ix[:max_ix]], \
                       node_no[ NX_ix[:max_ix]],\
                       mx[ MY_ix[:max_ix]], \
                       my[ MY_ix[:max_ix]], \
                       mxy[ MY_ix[:max_ix]], \
                       nx[ MY_ix[:max_ix]], \
                       ny[ MY_ix[:max_ix]], \
                       nxy[ MY_ix[:max_ix]], \
                       ###
                       m1[ NX_ix[:max_ix]], \
                       m2[ NX_ix[:max_ix]], \
                       alpha_M[ NX_ix[:max_ix]], \
                       mx_M[ NX_ix[:max_ix]],\
                       my_M[ NX_ix[:max_ix]], \
                       nx_M[ NX_ix[:max_ix]], \
                       ny_M[ NX_ix[:max_ix]], \
                       # NX
                       sig_n_NX[ NX_ix[:max_ix]], \
                       sig_m_NX[ NX_ix[:max_ix]], \
                       eta_n_NX[ NX_ix[:max_ix]], \
                       eta_m_NX[ NX_ix[:max_ix]], \
                       eta_tot_NX[ NX_ix[:max_ix]] )))
                                        
save_list_NX = save_arr_NX.tolist()

# specify the labels of the outputfile:
label_list_NX = ['row_no', \
              'elem_no', \
              'node_no', \
              'mx', \
              'my', \
              'mxy', \
              'nx' ,\
              'ny', \
              'nxy', \
              'n1', \
              'n2', \
              'alpha_N', \
              'mx_N', \
              'my_N', \
              'nx_N', \
              'ny_N', \
              'sig_n_NX', \
              'sig_m_NX', \
              'eta_n_NX', \
              'eta_m_NX', \
              'eta_tot_NX'
              ]

# save the array input to file 'report' applying formated strings
output_GdG_file = open('output_GdG','a')

caption_str_NX = str('GdG : Evaluation for case NX:')
label_str_NX = str()
dash_line_str = '\n'
for i in range(len(label_list_NX)):
    label_str_NX = label_str_NX + '%010s '
    dash_line_str = dash_line_str + '---------- '
    
header_str_NX = '\n \n \n' + caption_str_NX + dash_line_str + '\n' + label_str_NX +  dash_line_str + '\n'

output_GdG_file.write(header_str_NX %tuple(label_list_NX))

for line in save_list_NX:
    output_GdG_file.write('%10.0f %10.0f %10.0f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f \n' \
                     %(line[0], line[1], line[2], line[3], line[4], line[5], line[6], line[7], line[8], line[9], line[10], line[11], line[12], line[13], line[14], line[15], line[16], line[17], line[18], line[19], line[20] ) )

output_GdG_file.close()



# ------------------------------------------------------------
# NY - Save the output-file:
# ------------------------------------------------------------ 

max_ix = 4

save_arr_NY = transpose(vstack(( row_no[ NY_ix[:max_ix]], \
                       elem_no[ NY_ix[:max_ix]], \
                       node_no[ NY_ix[:max_ix]],\
                       mx[ NY_ix[:max_ix]], \
                       my[ NY_ix[:max_ix]], \
                       mxy[ NY_ix[:max_ix]], \
                       nx[ NY_ix[:max_ix]], \
                       ny[ NY_ix[:max_ix]], \
                       nxy[ NY_ix[:max_ix]], \
                       ###
                       m1[ NY_ix[:max_ix]], \
                       m2[ NY_ix[:max_ix]], \
                       alpha_M[ NY_ix[:max_ix]], \
                       mx_M[ NY_ix[:max_ix]],\
                       my_M[ NY_ix[:max_ix]], \
                       nx_M[ NY_ix[:max_ix]], \
                       ny_M[ NY_ix[:max_ix]], \
                       # NY
                       sig_n_NY[ NY_ix[:max_ix]], \
                       sig_m_NY[ NY_ix[:max_ix]], \
                       eta_n_NY[ NY_ix[:max_ix]], \
                       eta_m_NY[ NY_ix[:max_ix]], \
                       eta_tot_NY[ NY_ix[:max_ix]] )))
                                        
save_list_NY = save_arr_NY.tolist()

# specify the labels of the outputfile:
label_list_NY = ['row_no', \
              'elem_no', \
              'node_no', \
              'mx', \
              'my', \
              'mxy', \
              'nx' ,\
              'ny', \
              'nxy', \
              'n1', \
              'n2', \
              'alpha_N', \
              'mx_N', \
              'my_N', \
              'nx_N', \
              'ny_N', \
              'sig_n_NY', \
              'sig_m_NY', \
              'eta_n_NY', \
              'eta_m_NY', \
              'eta_tot_NY'
              ]

# save the array input to file 'report' applying formated strings
output_GdG_file = open('output_GdG','a')

caption_str_NY = str('GdG : Evaluation for case NY:')
label_str_NY = str()
dash_line_str = '\n'
for i in range(len(label_list_NY)):
    label_str_NY = label_str_NY + '%010s '
    dash_line_str = dash_line_str + '---------- '
    
header_str_NY = '\n \n \n' + caption_str_NY + dash_line_str + '\n' + label_str_NY +  dash_line_str + '\n'

output_GdG_file.write(header_str_NY %tuple(label_list_NY))

for line in save_list_NY:
    output_GdG_file.write('%10.0f %10.0f %10.0f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f \n' \
                     %(line[0], line[1], line[2], line[3], line[4], line[5], line[6], line[7], line[8], line[9], line[10], line[11], line[12], line[13], line[14], line[15], line[16], line[17], line[18], line[19], line[20] ) )

output_GdG_file.close()

