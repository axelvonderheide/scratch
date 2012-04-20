from numpy import *
from numpy.linalg import *
from yield_face import *
#from vtk_try import *
#from rm_to_vtk import *
#from yf_to_vtk import *

#STATE VARIABLES
epsilon_p_n = zeros((1,6))
delta_gamma = 0.
sigma_0 = 1.
f_diff1s = zeros((1,6))
f_diff1q = zeros((1,7)) 
f_diff2ss = zeros((6,6))  


#INPUT
epsilon_n = mat([0.,0.,0.,0.,0.,0.])
Deltas = (mat([0.,0.,0.,0.,0.,0.]), mat([1.,0.,0.,0.,0.,0.]),mat([1.0,0.,0.,0.,0.,0.]),mat([1.0,0.,0.,0.,0.,0.]), mat([-1.,0.,0.,0.,0.,0.]),mat([-1.0,0.,0.,0.,0.,0.]),mat([-1.0,0.,0.,0.,0.,0.]))
E = 1.0
nu= 0.2
K_bar = 0.  #parameter
H_bar = 0.  #parameter

h = mat([[K_bar,0.,0.,0.,0.,0.,0.],
         [0.,H_bar,0.,0.,0.,0.,0.],
         [0.,0.,H_bar,0.,0.,0.,0.],
         [0.,0.,0.,H_bar,0.,0.,0.],
         [0.,0.,0.,0.,H_bar,0.,0.],
         [0.,0.,0.,0.,0.,H_bar,0.],
         [0.,0.,0.,0.,0.,0.,H_bar]])   #the generalized plastic moduli

sigma_n = mat([0.,0.,0.,0.,0.,0.])
q_2 = mat([0.,0.,0.,0.,0.,0.])
q_1=0.0
max_iter=100

xi_trial = mat([0.,0.,0.,0.,0.,0.]).T
TOL = 1.e-7    # Tolerance

#Programm

es = 'e'
#es = 's'
#es = 0

#yf = J2()
yf = DruckerPrager()     
#yf = Gurson()
#yf = CamClay()

Cm = mat([[2*E/2/(1+nu)+E*nu/(1+nu)/(1-2*nu), E*nu/(1+nu)/(1-2*nu), E*nu/(1+nu)/(1-2*nu), 0, 0, 0],[E*nu/(1+nu)/(1-2*nu), 2*E/2/(1+nu)+E*nu/(1+nu)/(1-2*nu), E*nu/(1+nu)/(1-2*nu), 0, 0, 0],[E*nu/(1+nu)/(1-2*nu), E*nu/(1+nu)/(1-2*nu), 2*E/2/(1+nu)+E*nu/(1+nu)/(1-2*nu), 0, 0, 0],[0, 0, 0, E/2/(1+nu), 0, 0],[0, 0, 0, 0, E/2/(1+nu), 0],[0, 0, 0, 0, 0, E/2/(1+nu)]])


m = 0
for  d_epsilon in Deltas:
     m = m + 1
     
     epsilon_n = epsilon_n + d_epsilon
     sigma_n = Cm * (epsilon_n-epsilon_p_n).T
     xi_trial = sigma_n - q_2.T
  
     f_trial = yf.get_f_trial(xi_trial, q_1)
     print "f_trial", f_trial
     int_count = 0
     if f_trial < TOL:
          print "sigma_n", sigma_n
          #print "f_trial", f_trial
          #get_rm_elastic_vtk(m)
          #k_ela = sigma_n/epsilon_n
          #print "k_ela ", k_ela
          #get_diff1s = yf.get_diff1s(f_diff1s, epsilon_n, d_epsilon, epsilon_p_n, E, nu, q_1,q_2)
         
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
          while f_trial > TOL:
               if int_count == max_iter :#steps of interation less than the maximum
                    print "Maximal number of iteration reached"
                    break
               
    
               get_diff1s = yf.get_diff1s(f_diff1s, epsilon_n, d_epsilon, epsilon_p_n, E, nu, q_1,q_2)
               get_diff1q = yf.get_diff1q(f_diff1q, epsilon_n, d_epsilon, epsilon_p_n, E, nu, q_1,q_2)

               delta_gamma_2 = f_trial/(f_diff1s * Cm * f_diff1s.T + f_diff1q * h * f_diff1q.T )
               delta_gamma = delta_gamma + delta_gamma_2             
               epsilon_p_n = epsilon_p_n + delta_gamma_2 * f_diff1s
               
               q_1 = float(q_1 - delta_gamma_2 * K_bar * f_diff1q[0,0])
               q_2[0,0] = q_2 [0,0] - delta_gamma_2 * H_bar * f_diff1q[0,1]
               q_2[0,1] = q_2 [0,1] - delta_gamma_2 * H_bar * f_diff1q[0,2]
               q_2[0,2] = q_2 [0,2] - delta_gamma_2 * H_bar * f_diff1q[0,3]
               q_2[0,3] = q_2 [0,3] - delta_gamma_2 * H_bar * f_diff1q[0,4]
               q_2[0,4] = q_2 [0,4] - delta_gamma_2 * H_bar * f_diff1q[0,5]
               q_2[0,5] = q_2 [0,5] - delta_gamma_2 * H_bar * f_diff1q[0,6]
               
               sigma_n = Cm * (epsilon_n - epsilon_p_n).T
               xi_trial = sigma_n - q_2.T

               f_trial = yf.get_f_trial(xi_trial, q_1)
               print "f_trial_after", f_trial
              
               int_count = int_count + 1
               print "number of iterations", int_count
               
               if es == 'e': 
                    a = epsilon_n - epsilon_p_n # a is the engineering notation
                    b = diag(array([a[0,0],a[0,1],a[0,2]]),k=0)
                    b[1,2]=b[2,1]=a[0,3]
                    b[0,2]=b[2,0]=a[0,4]
                    b[0,1]=b[1,0]=a[0,5]

               elif es == 's':
                    a = sigma_n # a is the engineering notation
                    b = diag(array([a[0,0],a[1,0],a[2,0]]),k=0)
                    b[1,2]=b[2,1]=a[3,0]
                    b[0,2]=b[2,0]=a[4,0]
                    b[0,1]=b[1,0]=a[5,0]

               else:
                    print "FAILURE: Please choose es = 'e' or es = 's'"
               
               c = concatenate((c,eigvals(b)))
               
               #get_rm_vtk(E, nu,c,int_count,m,es)
          print "sigma_n ", sigma_n
          get_diff2ss = yf.get_diff2ss(f_diff2ss, epsilon_n, d_epsilon, epsilon_p_n, E, nu, q_2)
          
          theta = (Cm.I + delta_gamma[0,0] * f_diff2ss).I
          N = theta * f_diff1s.T/sqrt(f_diff1s*theta*f_diff1s.T)

          # stiffness
          #k_pla=theta-N*N.T
          #print "sigma_n", sigma_n  
          #print "k_pla ", k_pla
          #get_vtk(yf,E,nu,q_1,q_2,m,es)
     

