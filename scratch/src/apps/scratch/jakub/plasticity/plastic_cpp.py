from numpy import *
from numpy.linalg import *
from yield_face import *
#from yf_to_vtk import *
#from rm_to_vtk import *

# STATE VARIABLES
epsilon_p_n = zeros((1,6))
delta_gamma=0.
sigma_0 = 1.
                 
f_diff1s = zeros((1,6))
f_diff1q = zeros((1,7)) 
f_diff2ss = zeros((6,6))  
f_diff2sq = zeros((6,7))
f_diff2qq = zeros((7,7))

                 
# INPUT and set the initial values
epsilon_n = zeros((1,6))
sigma_n = zeros((1,6))

Deltas = (mat([0.,0.,0.,0.,0.,0.]), mat([1.,0.,0.,0.,0.,0.]),mat([1.,0.,0.,0.,0.,0.]),mat([1.,0.,0.,0.,0.,0.]), mat([-1.,0.,0.,0.,0.,0.]),mat([-1.,0.,0.,0.,0.,0.]),mat([-1.,0.,0.,0.,0.,0.]))

E  = 1.0
nu = 0.2

K_bar = 0.# parameter
H_bar = 0. # parameter

h = mat([[K_bar,0.,0.,0.,0.,0.,0.],
         [0.,H_bar,0.,0.,0.,0.,0.],
         [0.,0.,H_bar,0.,0.,0.,0.],
         [0.,0.,0.,H_bar,0.,0.,0.],
         [0.,0.,0.,0.,H_bar,0.,0.],
         [0.,0.,0.,0.,0.,H_bar,0.],
         [0.,0.,0.,0.,0.,0.,H_bar]])   # the generalized plastic moduli

max_iter = 100   
q_1 = 0.
q_2 = zeros((1,6))
TOL1 = 1.e-7   # Tolerance for the f_trial
#TOL2 = 1.e-5   # Tolerance for the R_n


# Choose one axis-system
es = 'e'    # set the 3 directions of strains als the three axises 
#es = 's'     # set the 3 directions of stresses als the three axises 
#es = 0      # lead to report of failure

# Choose one yield face
#yf = J2()       
yf = DruckerPrager()
#yf = Gurson()
#yf = Tresca()
#yf = Rankine()
#yf = CamClay()

 
# Programm  

Cm = mat([[2*E/2/(1+nu)+E*nu/(1+nu)/(1-2*nu), E*nu/(1+nu)/(1-2*nu), E*nu/(1+nu)/(1-2*nu), 0, 0, 0],[E*nu/(1+nu)/(1-2*nu), 2*E/2/(1+nu)+E*nu/(1+nu)/(1-2*nu), E*nu/(1+nu)/(1-2*nu), 0, 0, 0],[E*nu/(1+nu)/(1-2*nu), E*nu/(1+nu)/(1-2*nu), 2*E/2/(1+nu)+E*nu/(1+nu)/(1-2*nu), 0, 0, 0],[0, 0, 0, E/2/(1+nu), 0, 0],[0, 0, 0, 0, E/2/(1+nu), 0],[0, 0, 0, 0, 0, E/2/(1+nu)]])

bm = Cm.I
R_n = zeros((13,1))
#R_norm = 100*TOL2

d_e_q = zeros((1,13))
D_E_Q = zeros((1,13))
z1 = zeros((6,7)) 

q = 0
if H_bar:
      q = q + 6
if K_bar:
      q = q + 1

hn=[]  # hn is h's submatrix
if q == 1:
      hn= h[0:1,0:1]
elif q == 6:
      hn= h[1:7,1:7]
elif q == 7:
      hn= h

#print "hn", hn
Z1=z1[0:6,(7-q):7]      # Z1 is z1's submatrix
z2 = Z1.T

if q != 0:
      H = hn.I
      C_D = bmat('bm Z1; z2 H')  
else:
      H =[]
      C_D = bm   

f_diff1Q = f_diff1q[0:1,(7-q):7]
f_diff2SQ = f_diff2sq[0:6,(7-q):7]
f_diff2QQ = f_diff2qq[(7-q):7,(7-q):7]
R_N = R_n[0:(6+q),0:1]
#print "bm ", bm
m = 0
for d_epsilon in Deltas:
      m = m + 1
      
      epsilon_n = epsilon_n + d_epsilon
      sigma_n = Cm * (epsilon_n-epsilon_p_n).T
      xi_trial = sigma_n - q_2.T

      f_trial = yf.get_f_trial(xi_trial, q_1)
      print "f_trial", f_trial
      int_count = 0
      if f_trial < TOL1: 
            print "sigma_n", sigma_n
            #get_rm_elastic_vtk(m)
      else:  
            if es == 'e': 
                  a = epsilon_n - epsilon_p_n # a is the engineering notation
                  b = diag(array([a[0,0],a[0,1],a[0,2]]),k=0)  # b is tensorial notation
                  b[1,2]=b[2,1]=a[0,3]
                  b[0,2]=b[2,0]=a[0,4]
                  b[0,1]=b[1,0]=a[0,5]  # change the strain from the engineering notation with matrix form 6*1 into the tensorial notation with matrix form 3*3
              
            elif es == 's':
                  a = sigma_n # a is the engineering notation
                  b = diag(array([a[0,0],a[1,0],a[2,0]]),k=0)
                  b[1,2]=b[2,1]=a[3,0]
                  b[0,2]=b[2,0]=a[4,0]
                  b[0,1]=b[1,0]=a[5,0]
               
            else:
                  print "FAILURE: Please choose es = 'e' or es = 's'"
            
            c=eigvals(b) # get the eigenvector of matrix b
           
            #while f_trial > TOL1 or R_norm > TOL2:
            while f_trial > TOL1:
                  if int_count == max_iter :    # steps of interation less than the maximum
                        print "Maximal number of iteration reached"
                        break
                  
                  get_diff2ss = yf.get_diff2ss(f_diff2ss, epsilon_n, d_epsilon, epsilon_p_n, E, nu, q_2)
                  A_1 = bm + float(delta_gamma) * f_diff2ss
                  
                  get_diff2sq = yf.get_diff2sq(f_diff2sq, epsilon_n, d_epsilon, epsilon_p_n, E, nu, q_2)

                  get_diff2sq = yf.get_diff2sq(f_diff2sq, epsilon_n, d_epsilon, epsilon_p_n, E, nu, q_2)
                  f_diff2SQ = f_diff2sq[0:6,(7-q):7]
                  A_2 = float(delta_gamma) * f_diff2SQ
                  A_3 = float(delta_gamma) * (f_diff2SQ.T)

                  get_diff2qq = yf.get_diff2qq(f_diff2qq, epsilon_n, d_epsilon, epsilon_p_n, E, nu, q_2)
                  f_diff2QQ = f_diff2qq[(7-q):7,(7-q):7]
                  A_4 = H + float(delta_gamma) * f_diff2QQ  # calculate the four elements A_1,A_2,A_3 und A_4 of matrix A^(-1) seperately  
                               
                  b_n = bmat('A_1 A_2; A_3 A_4')  # combine the four elements into one matrix A^(-1)
                  A_n = b_n.I 
                 
                  get_diff1s = yf.get_diff1s(f_diff1s, epsilon_n, d_epsilon, epsilon_p_n, E, nu, q_1, q_2)
                  get_diff1q = yf.get_diff1q(f_diff1q, epsilon_n, d_epsilon, epsilon_p_n, E, nu, q_1, q_2)
                  f_diff1Q = f_diff1q[0:1,(7-q):7]
                  #print "f_diff1Q ", f_diff1Q
                  f_s_q = bmat('f_diff1s f_diff1Q')   # combine all the elements in f_diff1s and f_diff1q into one matrix f_s_q
                  #print "f_s_q ", f_s_q 

                  f_trial = yf.get_f_trial(xi_trial, q_1)
                  delta_gamma_2=(f_trial - f_s_q * A_n * R_N)/(f_s_q * A_n * f_s_q.T)  # calculate the increment to consistency parameter
                 
                  d_e_q = (C_D * A_n * (R_N + (delta_gamma_2 * f_s_q).T )).T  # obtain incremental plastic strains and internal variabes
                  D_E_Q[0,0:6]= d_e_q[0,0:6]
                  #print "d_e_q ", d_e_q
                  #print "C_D ", C_D
                  if q == 1:
                        D_E_Q[0,6]= d_e_q[0,6]
                  elif q == 6:
                        D_E_Q[0,-6:]= d_e_q[0,-6:]
                  elif q == 7:
                        D_E_Q[0,-7:]= d_e_q[0,-7:]

                  epsilon_p_n = epsilon_p_n + D_E_Q[0,0:6]
                  q_1 = -K_bar * D_E_Q[0,6] + q_1
                  q_2 = -H_bar * D_E_Q[0,-6:] + q_2
                 
                  sigma_n = Cm * (epsilon_n-epsilon_p_n).T
                  xi_trial = sigma_n - q_2.T
                  f_trial = yf.get_f_trial(xi_trial, q_1)
                  print "f_trial_after", f_trial
                
                  delta_gamma = delta_gamma + delta_gamma_2
                  
                  R_N = -d_e_q.T + (delta_gamma * f_s_q).T  # calculate the residual
                  R_norm = 0.
                  R_norm = sqrt(R_N.T * R_N) # obtain ||R_N||
                  int_count = int_count + 1
                  print "number of iterations", int_count

                  if es == 'e': 
                        a = epsilon_n - epsilon_p_n   # a is the engineering notation
                        b = diag(array([a[0,0],a[0,1],a[0,2]]),k=0)
                        b[1,2]=b[2,1]=a[0,3]
                        b[0,2]=b[2,0]=a[0,4]
                        b[0,1]=b[1,0]=a[0,5]

                  elif es == 's':
                        a = sigma_n   # a is the engineering notation
                        b = diag(array([a[0,0],a[1,0],a[2,0]]),k=0)
                        b[1,2]=b[2,1]=a[3,0]
                        b[0,2]=b[2,0]=a[4,0]
                        b[0,1]=b[1,0]=a[5,0]
                  
                  else:
                        print "FAILURE: Please choose es = 'e' or es = 's'"
               
                  #c = concatenate((c,eigvals(b)))
                  #get_rm_vtk(E, nu,c,int_count,m,es)
                  #print "c", c
                  #print "q_1 ",q_1,"q_2 ",q_2
              
            #get_diff2ss = yf.get_diff2ss(f_diff2ss, epsilon_n, d_epsilon, epsilon_p_n, E, nu, q_2)
           
            #theta = (Cm.I + delta_gamma[0,0] * f_diff2ss).I
            #N = theta * f_diff1s.T/sqrt(f_diff1s*theta*f_diff1s.T)

            # stiffness
            #k_pla=theta-N*N.T
            print "sigma_n", sigma_n  
            #print "k_pla ", k_pla
      #get_vtk(yf, E, nu, q_1, q_2, m, es)
