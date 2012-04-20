#from numpy import array, loadtxt, arange, sqrt, zeros, arctan, sin, cos, ones_like,\
#                zeros_like, vstack, savetxt, hstack, argsort, fromstring, transpose
from numpy import *
from math import pi
from string import split
from input_params import D_b, D_t, f_ctk, f_m, F_Rtex_l, F_Rtex_q


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
input_arr = loadtxt( 'adapted_data', usecols=(-8,-7,-6,-5,-4,-3), delimiter='\t' )

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

row_no = arange(elem_no.shape[0])


# ------------------------------------------------------------
# Parameters for the cracked state
# ------------------------------------------------------------
# derive the element thickness based on the element number
# (the division by 100 (=int) yields an integer)
D_elem = D_b - (elem_no/100-1)/11 * (D_b - D_t)

# (resultierende statische Nutzhoehe) 
d = 0.75 * D_elem

# (Abstand Schwereachse zur resultierende Bewehrungslage) 
# chose the same amount of reinforcement at the top as at the bottom 
# i.e. zs = zs1 = zs2
zs = d - D_elem/2.

# (Innerer Hebelarm) 
z = 0.9 * d


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
# Dimensioning
# ------------------------------------------------------------ 
# Index X: dimensioning in the transformed X-direction
# Index Y: dimensioning in the transformed Y-direction

# ------------------------------------------------------------
# case: MX
# ------------------------------------------------------------
# (Exzentrizitaet)
e_MX = mx_M / nx_M

# moment at the height of the resulting reinforcement layer:
m_Eds_MX = abs(mx_M) - zs * nx_M

# tensile force in the reinforcement for bending and compression
f_t_MX = m_Eds_MX / z + nx_M

# check if the two conditions are true:
zero_arr = zeros_like(nx)
cond_1 = nx_M > 0
cond_2 = e_MX < zs
bool_arr = cond_1 * cond_2
# in case of pure tension in the cross section:
f_t_MX[bool_arr] = nx_M * (zs + e_MX)/(zs + zs) 

# angel of deflection of the textile reinforcement (Umlenkwinkel)
beta = abs(alpha_M)

# resulting strength of the bi-directional textile considering the 
# deflection of the reinforcement in the loading direction:
f_Rtex_MX = F_Rtex_l * cos(pi/2 - beta)*(1-(pi/2 - beta))/(pi/2) + \
            F_Rtex_q * cos(beta)       *(1 -    beta    )/(pi/2) 

# necessary number of reinfocement layers
n_MX = f_t_MX / f_Rtex_MX

# sort the values depending on the maximum numer of reinfocement layers
MX_ix = argsort( n_MX )

# ------------------------------------------------------------
# case: MY
# ------------------------------------------------------------
# (Exzentrizitaet)
e_MY = my_M / ny_M

# moment at the height of the resulting reinforcement layer:
m_Eds_MY = abs(my_M) - zs * ny_M

# tensile force in the reinforcement for bending and compression
f_t_MY = m_Eds_MY / z + ny_M

# check if the two conditions are true:
zero_arr = zeros_like(nx)
cond_1 = nx_M > 0
cond_2 = e_MY < zs
bool_arr = cond_1 * cond_2
# in case of pure tension in the cross section:
f_t_MY[bool_arr] = ny_M * (zs + e_MY)/(zs + zs) 

# angel of deflection of the textile reinforcement (Umlenkwinkel)
beta = abs(alpha_M)

# resulting strength of the bi-directional textile considering the 
# deflection of the reinforcement in the loading direction:
f_Rtex_MY = F_Rtex_l * cos(pi/2 - beta)*(1-(pi/2 - beta))/(pi/2) + \
            F_Rtex_q * cos(beta)       *(1 -    beta    )/(pi/2) 

# necessary number of reinfocement layers
n_MY = f_t_MY / f_Rtex_MY

# sort the values depending on the maximum numer of reinfocement layers
MY_ix = argsort( n_MY )


# ------------------------------------------------------------
# case: NX
# ------------------------------------------------------------
# (Exzentrizitaet)
e_NX = mx_N / nx_N

# moment at the height of the resulting reinforcement layer:
m_Eds_NX = abs(mx_N) - zs * nx_N

# tensile force in the reinforcement for bending and compression
f_t_NX = m_Eds_NX / z + nx_N

# check if the two conditions are true:
zero_arr = zeros_like(nx)
cond_1 = nx_N > 0
cond_2 = e_NX < zs
bool_arr = cond_1 * cond_2
# in case of pure tension in the cross section:
f_t_NX[bool_arr] = nx_N * (zs + e_NX)/(zs + zs) 

# angel of deflection of the textile reinforcement (Umlenkwinkel)
beta = abs(alpha_N)

# resulting strength of the bi-directional textile considering the 
# deflection of the reinforcement in the loading direction:
f_Rtex_NX = F_Rtex_l * cos(pi/2 - beta)*(1-(pi/2 - beta))/(pi/2) + \
            F_Rtex_q * cos(beta)       *(1 -    beta    )/(pi/2) 

# necessary number of reinfocement layers
n_NX = f_t_NX / f_Rtex_NX

# sort the values depending on the maximum numer of reinfocement layers
NX_ix = argsort( n_NX )


# ------------------------------------------------------------
# case: NY
# ------------------------------------------------------------
# (Exzentrizitaet)
e_NY = ny_N / my_N

# moment at the height of the resulting reinforcement layer:
m_Eds_NY = abs(my_N) - zs * ny_N

# tensile force in the reinforcement for bending and compression
f_t_NY = m_Eds_NY / z + ny_N

# check if the two conditions are true:
zero_arr = zeros_like(nx)
cond_1 = nx_N > 0
cond_2 = e_NY < zs
bool_arr = cond_1 * cond_2
# in case of pure tension in the cross section:
f_t_NY[bool_arr] = ny_N * (zs + e_NY)/(zs + zs) 

# angel of deflection of the textile reinforcement (Umlenkwinkel)
beta = abs(alpha_N)

# resulting strength of the bi-directional textile considering the 
# deflection of the reinforcement in the loading direction:
f_Rtex_NY = F_Rtex_l * cos(pi/2 - beta)*(1-(pi/2 - beta))/(pi/2) + \
            F_Rtex_q * cos(beta)       *(1 -    beta    )/(pi/2) 

# necessary number of reinfocement layers
n_NY = f_t_NY / f_Rtex_NY

# sort the values depending on the maximum numer of reinfocement layers
NY_ix = argsort( n_NY )


# ------------------------------------------------------------
# Save the output-file:
# ------------------------------------------------------------ 

# ------------------------------------------------------------
# MX - Save the output-file:
# ------------------------------------------------------------ 

max_ix = 4

save_arr_MX = vstack(( \
                           #1
                           row_no[ MX_ix[:max_ix]],  \
                           #2
                           elem_no[ MX_ix[:max_ix]], \
                           #3
                           node_no[ MX_ix[:max_ix]], \
                           #4  
                           mx[ MX_ix[:max_ix]], \
                           #5
                           my[ MX_ix[:max_ix]], \
                           #6
                           mxy[ MX_ix[:max_ix]], \
                           #7
                           nx[ MX_ix[:max_ix]], \
                           #8
                           ny[ MX_ix[:max_ix]], \
                           #9
                           nxy[ MX_ix[:max_ix]], \
                           #10
                           m1[ MX_ix[:max_ix]], \
                           #12
                           m2[ MX_ix[:max_ix]], \
                           #13
                           alpha_M[ MX_ix[:max_ix]], \
                           #14  
                           mx_M[ MX_ix[:max_ix]],\
                           #15
                           my_M[ MX_ix[:max_ix]], \
                           #16
                           nx_M[ MX_ix[:max_ix]], \
                           #17
                           ny_M[ MX_ix[:max_ix]], \
                           #18
                           e_MX[ MX_ix[:max_ix]], \
                           #19
                           m_Eds_MX[ MX_ix[:max_ix]], \
                           #20
                           f_t_MX[ MX_ix[:max_ix]], \
                           #21
                           f_Rtex_MX[ MX_ix[:max_ix]], \
                           #22
                           n_MX[ MX_ix[:max_ix]] \
                           ))
                           

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
              'e_MX', \
              'm_Eds_MX', \
              'f_t_MX', \
              'f_Rtex_MX', \
              'n_MX'
              ]

# save the array input to file 'report' applying formated strings
output_GdT_file = open('output_GdT','w')

caption_str_MX = str('Evaluation for case MX:')
label_str_MX = str()
dash_line_str = '\n'

for i in range(len(label_list_MX)):
    label_str_MX = label_str_MX + '%010s '
    dash_line_str = dash_line_str + '---------- '
    
header_str_MX = caption_str_MX + dash_line_str + '\n' + label_str_MX +  dash_line_str + '\n'

output_GdT_file.write(header_str_MX %tuple(label_list_MX))

for line in save_list_MX:
    output_GdT_file.write('%10.0f %10.0f %10.0f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f \n' \
                     %(line[0], line[1], line[2], line[3], line[4], line[5], line[6], line[7], line[8], line[9], line[10], line[11], line[12], line[13], line[14], line[15], line[16], line[17], line[18], line[19], line[20] ) )

output_GdT_file.close()


## ------------------------------------------------------------
## MY - Save the output-file:
## ------------------------------------------------------------ 
#
#save_arr_MY = transpose(vstack(( \
#                           #1
#                           row_no[ MY_ix[:max_ix]],  \
#                           #2
#                           elem_no[ MY_ix[:max_ix]], \
#                           #3
#                           node_no[ MY_ix[:max_ix]], \
#                           #4  
#                           mx[ MY_ix[:max_ix]], \
#                           #5
#                           my[ MY_ix[:max_ix]], \
#                           #6
#                           mxy[ MY_ix[:max_ix]], \
#                           #7
#                           nx[ MY_ix[:max_ix]], \
#                           #8
#                           ny[ MY_ix[:max_ix]], \
#                           #9
#                           nxy[ MY_ix[:max_ix]], \
#                           #10
#                           m1[ MY_ix[:max_ix]], \
#                           #12
#                           m2[ MY_ix[:max_ix]], \
#                           #13
#                           alpha_M[ MY_ix[:max_ix]], \
#                           #14  
#                           mx_M[ MY_ix[:max_ix]],\
#                           #15
#                           my_M[ MY_ix[:max_ix]], \
#                           #16
#                           nx_M[ MY_ix[:max_ix]], \
#                           #17
#                           ny_M[ MY_ix[:max_ix]], \
#                           #18
#                           e_MY[ MY_ix[:max_ix]], \
#                           #19
#                           m_Eds_MY[ MY_ix[:max_ix]], \
#                           #20
#                           f_t_MY[ MY_ix[:max_ix]], \
#                           #21
#                           f_Rtex_MY[ MY_ix[:max_ix]], \
#                           #22
#                           n_MY[ MY_ix[:max_ix]] \
#                           )))
#                           
#
#save_list_MY = save_arr_MY.tolist()
#
## specify the labels of the outputfile:
#label_list_MY = ['row_no', \
#              'elem_no', \
#              'node_no', \
#              'mx', \
#              'my', \
#              'mxy', \
#              'nx' ,\
#              'ny', \
#              'nxy', \
#              'm1', \
#              'm2', \
#              'alpha_M', \
#              'mx_M', \
#              'my_M', \
#              'nx_M', \
#              'ny_M', \
#              'e_MY', \
#              'm_Eds_MY', \
#              'f_t_MY', \
#              'f_Rtex_MY', \
#              'n_MY'
#              ]
#
## save the array input to file 'output_GdT' applying formated strings
#output_GdT_file = open('output_GdT','a')
#
#caption_str = str('Evaluation for case MY:')
#label_str_MY = str()
#dash_line_str = '\n'
#for i in range(len(label_list_MY)):
#    label_str_MY = label_str_MY + '%010s '
#    dash_line_str = dash_line_str + '---------- '
#    
#header_str_MY = '\n \n \n ' + caption_str + dash_line_str + '\n' + label_str_MY  + dash_line_str + '\n'
#
#output_GdT_file.write(header_str_MY %tuple(label_list_MY))
#
#for line in save_list_MY:
#    output_GdT_file.write('%10.0f %10.0f %10.0f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f \n' \
#                     %(line[0], line[1], line[2], line[3], line[4], line[5], line[6], line[7], line[8], line[9], line[10], line[11], line[12], line[13], line[14], line[15], line[16], line[17], line[18], line[19], line[20] ) )
#
#output_GdT_file.close()
#
#
#
## ------------------------------------------------------------
## NX - Save the output-file:
## ------------------------------------------------------------ 
#
#save_arr_NX = transpose(vstack(( \
#                           #1
#                           row_no[ NX_ix[:max_ix]],  \
#                           #2
#                           elem_no[ NX_ix[:max_ix]], \
#                           #3
#                           node_no[ NX_ix[:max_ix]], \
#                           #4  
#                           mx[ NX_ix[:max_ix]], \
#                           #5
#                           my[ NX_ix[:max_ix]], \
#                           #6
#                           mxy[ NX_ix[:max_ix]], \
#                           #7
#                           nx[ NX_ix[:max_ix]], \
#                           #8
#                           ny[ NX_ix[:max_ix]], \
#                           #9
#                           nxy[ NX_ix[:max_ix]], \
#                           #10
#                           n1[ NX_ix[:max_ix]], \
#                           #12
#                           n2[ NX_ix[:max_ix]], \
#                           #13
#                           alpha_N[ NX_ix[:max_ix]], \
#                           #14  
#                           mx_N[ NX_ix[:max_ix]],\
#                           #15
#                           my_N[ NX_ix[:max_ix]], \
#                           #16
#                           nx_N[ NX_ix[:max_ix]], \
#                           #17
#                           ny_N[ NX_ix[:max_ix]], \
#                           #18
#                           e_NX[ NX_ix[:max_ix]], \
#                           #19
#                           m_Eds_NX[ NX_ix[:max_ix]], \
#                           #20
#                           f_t_NX[ NX_ix[:max_ix]], \
#                           #21
#                           f_Rtex_NX[ NX_ix[:max_ix]], \
#                           #22
#                           n_NX[ NX_ix[:max_ix]] \
#                           )))
#                           
#
#save_list_NX = save_arr_NX.tolist()
#
## specify the labels of the outputfile:
#label_list_NX = ['row_no', \
#              'elem_no', \
#              'node_no', \
#              'mx', \
#              'my', \
#              'mxy', \
#              'nx' ,\
#              'ny', \
#              'nxy', \
#              'n1', \
#              'n2', \
#              'alpha_N', \
#              'mx_N', \
#              'my_N', \
#              'nx_N', \
#              'ny_N', \
#              'e_NX', \
#              'm_Eds_NX', \
#              'f_t_NX', \
#              'f_Rtex_NX', \
#              'n_NX'
#              ]
#
## save the array input to file 'output_GdT' applying formated strings
#output_GdT_file = open('output_GdT','a')
#
#caption_str = str('Evaluation for case NX: ')
#label_str_NX = str()
#dash_line_str = '\n'
#for i in range(len(label_list_NX)):
#    label_str_NX = label_str_NX + '%010s '
#    dash_line_str = dash_line_str + '---------- '
#    
#header_str_NX = '\n \n \n ' + caption_str + dash_line_str + '\n' + label_str_NX + dash_line_str + '\n'
#
#output_GdT_file.write(header_str_NX %tuple(label_list_NX))
#
#for line in save_list_NX:
#    output_GdT_file.write('%10.0f %10.0f %10.0f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f \n' \
#                     %(line[0], line[1], line[2], line[3], line[4], line[5], line[6], line[7], line[8], line[9], line[10], line[11], line[12], line[13], line[14], line[15], line[16], line[17], line[18], line[19], line[20] ) )
#
#output_GdT_file.close()
#
#
## ------------------------------------------------------------
## NY - Save the output-file:
## ------------------------------------------------------------ 
#
#save_arr_NY = transpose(vstack(( \
#                           #1
#                           row_no[ NY_ix[:max_ix]],  \
#                           #2
#                           elem_no[ NY_ix[:max_ix]], \
#                           #3
#                           node_no[ NY_ix[:max_ix]], \
#                           #4  
#                           mx[ NY_ix[:max_ix]], \
#                           #5
#                           my[ NY_ix[:max_ix]], \
#                           #6
#                           mxy[ NY_ix[:max_ix]], \
#                           #7
#                           nx[ NY_ix[:max_ix]], \
#                           #8
#                           ny[ NY_ix[:max_ix]], \
#                           #9
#                           nxy[ NY_ix[:max_ix]], \
#                           #10
#                           n1[ NY_ix[:max_ix]], \
#                           #12
#                           n2[ NY_ix[:max_ix]], \
#                           #13
#                           alpha_N[ NY_ix[:max_ix]], \
#                           #14  
#                           mx_N[ NY_ix[:max_ix]],\
#                           #15
#                           my_N[ NY_ix[:max_ix]], \
#                           #16
#                           nx_N[ NY_ix[:max_ix]], \
#                           #17
#                           ny_N[ NY_ix[:max_ix]], \
#                           #18
#                           e_NY[ NY_ix[:max_ix]], \
#                           #19
#                           m_Eds_NY[ NY_ix[:max_ix]], \
#                           #20
#                           f_t_NY[ NY_ix[:max_ix]], \
#                           #21
#                           f_Rtex_NY[ NY_ix[:max_ix]], \
#                           #22
#                           n_NY[ NY_ix[:max_ix]] \
#                           )))
#                           
#
#save_list_NY = save_arr_NY.tolist()
#
## specify the labels of the outputfile:
#label_list_NY = ['row_no', \
#              'elem_no', \
#              'node_no', \
#              'mx', \
#              'my', \
#              'mxy', \
#              'nx' ,\
#              'ny', \
#              'nxy', \
#              'm1', \
#              'm2', \
#              'alpha_M', \
#              'mx_M', \
#              'my_M', \
#              'nx_M', \
#              'ny_M', \
#              'e_NY', \
#              'm_Eds_NY', \
#              'f_t_NY', \
#              'f_Rtex_NY', \
#              'n_NY'
#              ]
#
## save the array input to file 'output_GdT' applying formated strings
#output_GdT_file = open('output_GdT','a')
#
#caption_str = str('Evaluation for case NY:')
#label_str_NY = str()
#dash_line_str = '\n'
#for i in range(len(label_list_NY)):
#    label_str_NY = label_str_NY + '%010s '
#    dash_line_str = dash_line_str + '---------- '
#    
#header_str_NY = '\n \n \n ' + caption_str + dash_line_str + '\n' + label_str_NY + dash_line_str + '\n'
#
#output_GdT_file.write(header_str_NY %tuple(label_list_NY))
#
#for line in save_list_NY:
#    output_GdT_file.write('%10.0f %10.0f %10.0f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f \n' \
#                     %(line[0], line[1], line[2], line[3], line[4], line[5], line[6], line[7], line[8], line[9], line[10], line[11], line[12], line[13], line[14], line[15], line[16], line[17], line[18], line[19], line[20] ) )
#
#output_GdT_file.close()

# save array to file 'output'
output_file = open('output_02','a')
savetxt('output_02',save_arr_MX.T, fmt = '%10.0f')
# ------------------
